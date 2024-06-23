

namespace ffd::combination::unit_test{


  bool combination(){
    bool IsOk = true;



    constexpr int k = 3;

    std::vector<std::array<int, k>> values_vec;
    
    for( auto values: ffd::user_space::
	   Combination<k>(ffd::vector_range::Range(3)) ){
      values_vec.push_back(values);
      // for( int j: ffd::vector_range::Range(k) ){
      // 	std::cerr<<values[j]<<" ";
      // }
      // std::cerr<<std::endl;
    }

    for(  int j: ffd::vector_range::Range( std::size(values_vec) )  ){
      switch(j){
      case 0:
	IsOk = IsOk && values_vec[j] == std::array{0, 0, 0};
	break;
      case 1:
	IsOk = IsOk && values_vec[j] == std::array{0, 0, 1};
	break;
      case 2:
	IsOk = IsOk && values_vec[j] == std::array{0, 1, 1};
	break;
      case 3:
	IsOk = IsOk && values_vec[j] == std::array{1, 1, 1};
	break;
      case 4:
	IsOk = IsOk && values_vec[j] == std::array{0, 0, 2};
	break;
      case 5:
	IsOk = IsOk && values_vec[j] == std::array{0, 1, 2};
	break;
      case 6:
	IsOk = IsOk && values_vec[j] == std::array{1, 1, 2};
	break;
      case 7:
	IsOk = IsOk && values_vec[j] == std::array{0, 2, 2};
	break;
      case 8:
	IsOk = IsOk && values_vec[j] == std::array{1, 2, 2};
	break;
      case 9:
	IsOk = IsOk && values_vec[j] == std::array{2, 2, 2};
	break;
      case 10:
	IsOk = false;
	break;
      }
    }

    return IsOk;
  }

}//namespace
