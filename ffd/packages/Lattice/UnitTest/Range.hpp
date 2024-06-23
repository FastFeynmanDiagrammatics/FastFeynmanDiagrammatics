namespace ffd::lattice::unit_test{

  void range(){

    int const Lx = 5, Ly = 7, Lz = 3;

    
    HypercubicLattice<3> H({Lx, Ly, Lz});

    {
      [[maybe_unused]] auto [vector_coordinates,
			     vector_int] = Range(H);

    
      assert(( size( vector_coordinates ) == size( vector_int ) ));
      assert(( size( vector_coordinates ) == Lx*Ly*Lz ));

      // for( auto coordinates_int: vector_int ){
      //   for( auto component: coordinates_int ){
      // 	std::cerr<<component<<" ";
      //   }
      //   std::cerr<<std::endl;
      // }

      using ffd::vector_range::Range;
      for( auto j: Range( size(vector_int) ) ){
	auto X =       vector_coordinates[j];
	auto X_array = vector_int[j];


	int x0 = ffd::get<0>(X)();
	x0 = (x0 >= 0 ? x0 : x0+Lx);
	assert(( x0 == X_array[0] ));


	int x1 = ffd::get<1>(X)();
	// std::cerr<<x1<<" ";
	x1 = (x1 >= 0 ? x1 : x1+Ly);
	// std::cerr<<x1<<" "<<X_array[1]<<" "<<std::endl;
	assert(( x1 == X_array[1] ));


	int x2 = ffd::get<2>(X)();
	x2 = (x2 >= 0 ? x2 : x2+Lz);
	assert(( x2 == X_array[2] ));
      }
    }




    
    KagomeLattice K({Lx, Ly});

    {
      [[maybe_unused]] auto [vector_coordinates,
			     vector_int] = Range(K);

    
      assert(( size( vector_coordinates ) == size( vector_int ) ));
      assert(( size( vector_coordinates ) == Lx*Ly*3 ));

      // for( auto coordinates_int: vector_int ){
      //   for( auto component: coordinates_int ){
      // 	std::cerr<<component<<" ";
      //   }
      //   std::cerr<<std::endl;
      // }
      
      
      using ffd::vector_range::Range;
      for( auto j: Range( size(vector_int) ) ){
	auto X =       vector_coordinates[j];
	auto X_array = vector_int[j];


	int x0 = ffd::get<0>(X)();
	x0 = (x0 >= 0 ? x0 : x0+Lx);
	assert(( x0 == X_array[0] ));


	int x1 = ffd::get<1>(X)();
	// std::cerr<<x1<<" ";
	x1 = (x1 >= 0 ? x1 : x1+Ly);
	// std::cerr<<x1<<" "<<X_array[1]<<" "<<std::endl;
	assert(( x1 == X_array[1] ));


	int x2 = ffd::get<2>(X)();
	assert(( x2 == X_array[2] ));
      }
    }    

    
  }

}//namespace
