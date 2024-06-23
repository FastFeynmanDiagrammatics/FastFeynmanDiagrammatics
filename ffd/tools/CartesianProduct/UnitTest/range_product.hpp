

namespace ffd::cartesian_product::unit_test{

  void range_product(){

    using namespace ffd::user_space;
    using ffd::vector_range::Range;
    
    constexpr int d = 3;
    
    std::array<int, d> L = {1, 3, 2};

    std::vector<std::array<int, d>> history_1;
    for( auto [x, y, z]: CartesianProductRange(L) ){
      // std::cerr<<x<<" "<<y<<" "<<z<<std::endl;
      history_1.push_back({x, y, z});
    }

    

    for(  int j: Range( size(history_1) )  ){
      switch(j){
      case 0:
	assert(( history_1[j] == std::array{0, 0, 0} ));
	break;
      case 1:
	assert(( history_1[j] == std::array{0, 0, 1} ));
	break;
      case 2:
	assert(( history_1[j] == std::array{0, 1, 0} ));
	break;
      case 3:
	assert(( history_1[j] == std::array{0, 1, 1} ));
	break;
      case 4:
	assert(( history_1[j] == std::array{0, 2, 0} ));
	break;
      case 5:
	assert(( history_1[j] == std::array{0, 2, 1} ));
	break;
      case 6:
	assert(false);
      };
    }
    assert(( size(history_1) == 6 ));
    


    // std::cerr<<"here"<<std::endl;
    std::vector<std::array<int, 2>> history_2;
    for( auto [x, y]: CartesianProductRange<2>(1, 3) ){
      // std::cerr<<x<<" "<<y<<std::endl;
      history_2.push_back({x, y});
    }


    // for( auto v: CartesianProductRange<2>(1, 3) ){
    //   std::cerr<<v[0]<<" "<<v[1]<<std::endl;
    //   history_2.push_back(v);
    // }


    assert(size(history_2) == 4);
    for( int j: Range( size(history_2) ) ){
      switch(j){
      case 0:
	assert(( history_2[j] == std::array<int, 2>{1, 1} ));
	break;
      case 1:
	assert(( history_2[j] == std::array<int, 2>{1, 2} ));
	break;
      case 2:
	assert(( history_2[j] == std::array<int, 2>{2, 1} ));
	break;
      case 3:
	assert(( history_2[j] == std::array{2, 2} ));
	break;
      case 4:
	assert(false);
      }
    }

    
  }


}//namespace
