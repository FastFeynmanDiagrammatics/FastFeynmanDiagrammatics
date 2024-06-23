

namespace ffd::cartesian_product::unit_test{

  void try_to_assert_example(){
    using namespace ffd::user_space;
    
    std::set<int> S1{{3, -2, 1}};
    std::set<int> S2{{4, 5}};
    std::set<int> S3 = {1};
    std::set<int> S4 = {2};
    std::array<std::set<int>, 4> array_sets{{S1, S2, S3, S4}};
    

    std::vector<std::array<int, 4>> rec_coordinates;

    {
      auto UU = ffd::user_space::CartesianProduct<4, std::set, int>(array_sets);
      ++UU;
    }
    int test1;

    // std::cerr<<"Here";
    for( auto [x, y, z, w]: ffd::user_space::CartesianProduct<4, std::set, int>(array_sets) ){
      // std::cerr<<x<<" "<<y<<" "<<z<<" "<<w<<std::endl;
      rec_coordinates.push_back({x, y, z, w});
    }

    

    
    for(  int j: ffd::vector_range::Range( size(rec_coordinates) )  ){
      switch(j){
      case 0:
	assert(( rec_coordinates[j] == std::array{-2, 4, 1, 2} ));
	break;
      case 1:
	assert(( rec_coordinates[j] == std::array{-2, 5, 1, 2} ));
	break;
      case 2:
	assert(( rec_coordinates[j] == std::array{1, 4, 1, 2} ));
	break;
      case 3:
	assert(( rec_coordinates[j] == std::array{1, 5, 1, 2} ));
	break;
      case 4:
	assert(( rec_coordinates[j] == std::array{3, 4, 1, 2} ));
	break;
      case 5:
	assert(( rec_coordinates[j] == std::array{3, 5, 1, 2} ));
	break;
      case 6:
	assert((false));
      };
    }


  }

}//namespace
