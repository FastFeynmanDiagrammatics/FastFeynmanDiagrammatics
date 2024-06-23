

namespace ffd::bubbles{

  // template<typename ImaginaryTime_t,
  // 	   typename Lattice_t>

  // void
  
  // Bubbles<ImaginaryTime_t,
  // 	  Lattice_t>::
  // CookLobster(Real precision){
  //   using ffd::vector_range::Range;
  //   if( precision < 0 ){
  //     precision = AbsolutePrecision;
  //   }
    

  //   auto Bravais_sites = ffd::user_space::
  //     CartesianProductRange(L_array);
  //   auto Cell_sites_x2 = ffd::user_space::
  //     CartesianProductRange<2*2>(cell_size);

    
  //   Lobster = std::vector<Lobster_bravais_t>(L_product*L_product);

    
  //   for( auto r2: ffd::user_space::
  // 	   CartesianPower<2>(Bravais_sites) ){
  //     std::array<std::array<int, d>, 2> r_array;
  //     for( int spin: {0, 1} ){
  // 	for( int u: Range(d) ){
  // 	  r_array[spin][u] = r2[u+d*spin];
  // 	}
  //     }
      
  //     std::array<int, 2*d> L_array_x2;
  //     for( int j: Range(2*d) ){
  // 	L_array_x2[j] = L_array[j%d];
  //     }

  //     int r2_i = create_index(r2, L_array_x2);
      
      
  //     for( auto n2: Cell_sites_x2){
  // 	std::array<std::array<int, 2>, 2> n_array =
  // 	  DistributeArray<cell_size>(;

	
  // 	std::array<int, 2*2> cell_size_x2_array;
  // 	cell_size_x2_array.fill(cell_size);
  // 	int n2_i = create_index(n2, cell_size_x2_array);

	
  // 	Lobster[r2_i][n2_i] = chebyshev_2d_t();

  // 	for( auto y: Bravais_sites){
  // 	  int y_i = create_index(y, L_array);
	  

  // 	  std::array<int, 2> r_y;
  // 	  for( int spin: {0, 1} ){
  // 	    auto r_y[spin] = substraction_indices<d>(r_array[spin],
  // 						     y,
  // 						     L_array);
  // 	  }
	  
	  
  // 	  for( auto m: Cell_sites){
  // 	    int m_i = create_index(m, std::array{cell_size, cell_size});


	    



  // 	  }
  // 	}
  //     }
      
  //   }
    

    
  // }


}//namespace
