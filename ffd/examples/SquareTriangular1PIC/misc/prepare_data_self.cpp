#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <array>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>

using Real = long double;
using ulong = unsigned long;

#include <ffd/tools/UncertainNumber.hpp>
#include <ffd/tools/WeightedMean.hpp>
#include <filesystem>

using wmean_t = ffd::user_space::WeightedMean<Real>;
using val_err_t = ffd::user_space::ENumber<Real>;

using namespace std; using namespace ffd::user_space;
namespace fs = std::filesystem;
constexpr Real epsilon = 1e-200;
constexpr int max_order = 32;

int main(int argc, char* argv[]){
  bool no_hold = true, verbose = false;
  std::map<std::string, std::string> file_extension;
  file_extension["raw standard series"] = ".s";
  //// format raw standard series
  // -1 normalization normalization_err
  // 0 normalized_order_0 normalized_order_0_err
  // 1 unnormalized_order_1 unnormalized_order_1_err
  // ...
  file_extension["standard series"] = ".series";
  //// format standard series
  // 0 order_0 order_0_err
  // 1 order_1 order_1_err
  // ...
  file_extension["raw double series"] = ".ds";
  //// format raw double series
  // -1 -1 normalization normalization_err
  // 0 0 normalized_order_{0, 0} normalized_order_{0,0}_err
  // ...
  // j k unnormalized_order_{j, k} unnormalized_order_{j,k}_err
  // ...
  file_extension["double series"] = ".dseries";
  //// format double series
  // 0 0 order_{0, 0} order_{0,0}_err
  // ...
  // j k order_{j, k} order_{j,k}_err
  // ...
  std::vector<std::string> paths;
  for ( ulong j = 1; j < argc; ++j ) {
    if ( !no_hold ) 
      ;
    else {
      if ( argv[j] == "-v" ) {
	verbose = true;
      } else {
	paths.push_back(argv[j]);
      }
    }
  } // for j in range(1, argc)

  
  auto compute_mean_dir =
    [&file_extension, verbose]
    (fs::path dir_path) mutable
    {
      std::string path_str = dir_path;
      if ( path_str.back() == '/' ) {
	path_str = path_str.substr(0, size(path_str)-1);
      }
      //      path_str.erase(std::remove(path_str.begin(), path_str.end(), '/'), path_str.end());
      if ( !fs::exists(dir_path) ) {
	std::cerr << "!!!!!!@@@@@@@#####$$$$$$$ " << dir_path << " does not exists!!!!!\n";
	return;
      }
      auto const current_path = fs::current_path();
      fs::current_path(dir_path);
      std::map<std::pair<int, int>, std::array<wmean_t, max_order>> mp;
      std::map<std::tuple<int, int, int>, wmean_t> mp_d; // TODO
      for ( auto& data_file: fs::directory_iterator(".") ) {
	if ( fs::is_regular_file(data_file.path()) &&
	     data_file.path().extension() == file_extension["raw standard series"] ) {
	  //	  std::cerr << "data_file_path = " << data_file.path() << "\n";
	  std::string fname_stem = data_file.path().stem();
	  int pid =
	    [fname_stem]
	    {
	      std::istringstream iss(fname_stem);
	      std::string pid_string;
	      std::getline(iss, pid_string, '_');
	      return std::stoi(pid_string);
	    }();
	  //	  std::cerr << "pid = " << pid;
	  auto const [val, order] =
	    [data_file]
	    {
	      std::ifstream df (data_file.path());
	      std::array<wmean_t, max_order> val;
	      Real value, error;
	      int order, order_temp;
	      while ( df >> order_temp >> value >> error ) {
		order = order_temp;
		val[order+1] = wmean_t(value, error);
	      }
	      return std::make_pair(val, order);
	    }();
	  //	  std::cerr << ", order = " << order << "\n";
	  auto const key = std::pair<int, int>(pid, order);
	  for ( ulong j = 0; j < order+2; ++j ) {
	    mp[key][j] &= val[j];
	  } // for j in range(0, order+2)
	} // if raw standard series
	if ( fs::is_regular_file(data_file.path()) &&
	     data_file.path().extension() == file_extension["raw double series"] ) {
	  std::string fname_stem = data_file.path().stem();
	  int pid =
	    [fname_stem]
	    {
	      std::istringstream iss(fname_stem);
	      std::string pid_string;
	      std::getline(iss, pid_string, '_');
	      return std::stoi(pid_string);
	    }();
	  //	  std::cerr << "pid = " << pid;
	  {
	    std::ifstream df (data_file.path());
	    Real value, error;
	    int order0, order1;
	    while ( df >> order0 >> order1 >> value >> error ) {
	      mp_d[std::tuple<int, int, int>{pid, order0, order1}] &= wmean_t(value, error);
	    }
	  }
	} // if raw double series
      } // for data_file in directory
      fs::current_path(current_path);
      if ( !mp.empty() ) {
	int order_max = 0;
	std::array<wmean_t, max_order> merged;
	for ( auto [key, val] : mp ) {
	  auto [pid, order] = key;
	  order_max = std::max(order, order_max);
	  merged[0] = val[1];
	  for ( ulong j = 1; j < order+1; ++j ) {
	    auto const val_j = val_err_t(val[j+1].eval()) / val_err_t(val[0].eval());
	    merged[j] &= wmean_t(val_j.eval());
	  } // for j in range(0, order+1)
	}
	if ( verbose )
	  std::cerr << dir_path << " merged\n";
	std::ofstream out (path_str+file_extension["standard series"]);
	for ( ulong j = 0; j < order_max+1; ++j ) {
	  auto [val, err] = merged[j].eval();
	  out << j << " " << val << " " << err << "\n";
	} // for j in range(0, order_max+1)
      } else if ( !mp_d.empty() ) {
	std::map<std::pair<int, int>, wmean_t> merged;
	for ( auto [key, val] : mp_d ) {
	  auto [pid, order0, order1] = key;
	  if ( order0 != -1 && order1 != -1 ) {
	    auto val_temp = val_err_t(val.eval());
	    if (order0 != 0 || order1 != 0 )
	      val_temp = val_temp / val_err_t(mp_d[std::make_tuple(pid, -1, -1)].eval());
	    merged[std::pair<int, int>{order0, order1}] &= wmean_t(val_temp.eval());

	  }
	}
	if ( verbose )
	  std::cerr << dir_path << " merged\n";
	std::ofstream out (path_str+file_extension["double series"]);
	for ( auto [key, val] : merged ) {
	  auto const [order0, order1] = key;
	  auto const [value, err] = val.eval();
	  out << order0 << " " << order1 << " " << value << " " << err << "\n";
	} // for key, val in merged
      }
    };

  auto nested_mean_dir =
    [&compute_mean_dir, verbose]
    (fs::path path) mutable -> void
    {
      //      std::cerr << "path = " << path << "\n";
      auto impl =
	[&compute_mean_dir, verbose]
	(fs::path path, auto& mean_dir) mutable -> void
	{
	  if ( verbose )
	    std::cerr << "path = " << path << "\n";
	  if ( fs::is_directory(path) && path != "./lambda_norm" && path != "./parameters" ) {
	    compute_mean_dir(path);
	    fs::current_path(path);
	    for ( auto& path_nested :
		    fs::directory_iterator(".") ) {
	      mean_dir(path_nested, mean_dir);
	    }
	    fs::current_path("..");
	  }
	};
      
      return impl(path, impl);
    };

  
  for ( auto dir_path : paths ) {
    nested_mean_dir(dir_path);
  }

  
} // main
  
