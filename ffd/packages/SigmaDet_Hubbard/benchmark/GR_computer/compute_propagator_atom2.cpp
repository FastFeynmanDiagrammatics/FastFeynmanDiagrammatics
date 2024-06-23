#include<ffd/std.hpp>
#include<ffd/core.hpp>
#include<ffd/tools.hpp>
#include<ffd/packages.hpp>
#include"parameters.hpp"
#include"G0.hpp"

using namespace ffd::user_space;



int main(){
  std::cerr<<std::setprecision(15)<<"nR = "<<-2*GR(Beta)<<std::endl;
  std::cerr<<"n0 = "<<-2*G0(Beta)<<std::endl;
  std::ostringstream out_data_str, file_sigma_name, file_G_name, common_part, benchmark_signature_name;
  // int Beta_integer = int(Beta+1e-10);
  // int U_integer = int(U*1+1e-10);
  out_data_str<<"out_data/";
  common_part<<"_B"<<Beta;
  common_part<<"_U"<<U;
  common_part<<"_1x1_";
  common_part<<"Nc"<<N_Chebyshev_DMFT;
  common_part<<"_name";
  common_part<<"Atom";
  common_part<<"_m"<<std::setprecision(15)<<mu;
  common_part<<"_A"<<std::setprecision(15)<<alpha_shift;
  common_part<<"_D"<<std::setprecision(2)<<Delta;
  common_part<<"_B"<<std::setprecision(3)<<Beta_Sigma;
  common_part<<".dat";
  file_sigma_name<<out_data_str.str();
  file_G_name<<out_data_str.str();
  file_sigma_name<<"SigmaR";
  file_G_name    <<"GR____";
  file_sigma_name<<common_part.str();
  file_G_name<<common_part.str();
  
  
  std::ofstream sigma_dmft_file( file_sigma_name.str() );
  std::ofstream g_dmft_file( file_G_name.str() );

  
  for( int j = n_frequencies_DMFT-1; j >= -n_frequencies_DMFT; --j ){
    sigma_dmft_file << std::setprecision(15)<<std::scientific<<matsu_freq(j)<<"\t";
    auto sigma_R = Sigma_R_omega(j);
    sigma_dmft_file << std::real(sigma_R)<<"\t";
    sigma_dmft_file << std::imag(sigma_R)<<"\n";
  }

  
  for(int j=0; j<N_Chebyshev_DMFT; ++j){
    using std::cos;
    using ffd::core_math::Pi;
    Real const tau_j = .5*Beta*(1.-cos((.5+j)*Pi/N_Chebyshev_DMFT));
    g_dmft_file<<std::setprecision(15)<<std::scientific<<
      "0\t"<<"0\t"<<j<<"\t"<<tau_j<<"\t"<<GR(tau_j)<<"\n";
  }

  
}

