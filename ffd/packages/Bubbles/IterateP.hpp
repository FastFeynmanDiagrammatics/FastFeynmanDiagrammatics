namespace ffd::bubbles{

  template<typename ImaginaryTime_t,
	   typename Lattice_t>

  typename Bubbles<ImaginaryTime_t,
		   Lattice_t>::error_t
  
  Bubbles<ImaginaryTime_t,
	  Lattice_t>::
  IterateP(){
    Real difference_between_iterations = 0.;

    using ffd::vector_range::Range;
    using std::size, std::abs;


    Gamma_t RPA_Bubble(L_product);
    auto Bravais_sites = ffd::user_space::
      CartesianProductRange(L_array);
    auto Cell_sites = ffd::user_space::
      CartesianProductRange<2>(cell_size);
      
    
    
    for( auto r: Bravais_sites ){
      int r_i = create_index(r, L_array);

      for( auto n: Cell_sites ){
	int n_i = create_index(n, std::array{cell_size, cell_size});

	
	auto RPA_Bubble_func = [r_i, n_i, this](Real tau)->Real{
				 if( RPA == RPA_scheme::ph ){
				   return +G[0][r_i][n_i](tau)*
				     G[1][r_i][n_i](beta-tau);
				 }else if (RPA == RPA_scheme::pp){
				   return -G[0][r_i][n_i](tau)*
				     G[1][r_i][n_i](tau);
				 }else{
				   assert((false));
				   return 0;
				 }
			       };
      

	RPA_Bubble[r_i][n_i] =
	  chebyshev_t(RPA_Bubble_func,
		      {0, beta},
		      AbsolutePrecision/U/U);
      
      }
    }

    
    Gamma_t Int_Bubble_P(L_product);
    for( auto r: Bravais_sites ){
      for( auto s: Bravais_sites ){
	auto r_s = substraction_indices<d>(r,
					   s,
					   L_array);
	int r_i = create_index(r, L_array);
	int s_i = create_index(s, L_array);
	int r_s_i = create_index(r_s, L_array);

	for( auto n: Cell_sites ){
	  for( auto m: Cell_sites ){
	    std::array<int, 2> cell_sizes{cell_size, cell_size};
	    static_assert( std::is_same_v<decltype(n), std::array<int, 2>> );
	    static_assert( std::is_same_v<decltype(m), std::array<int, 2>> );
	    static_assert( std::is_same_v<decltype(cell_sizes), std::array<int, 2>> );
	    auto n_m = substraction_indices<2>(n,
					       m,
					       cell_sizes);
	    int n_i = create_index(n, cell_sizes);
	    int m_i = create_index(m, cell_sizes);
	    int n_m_i = create_index(n_m, cell_sizes);
	  
	    Int_Bubble_P[r_i][n_i] += ffd::imaginary_time_convolution::
	      Convolute1And2FromZeroTo3<Real,
					chebyshev_t>(RPA_Bubble[s_i][m_i],
						     P[r_s_i][n_m_i],
						     beta,
						     false,//bosonic
						     AbsolutePrecision/U/U);
// ffd::convolution::
// 	      Convolution(RPA_Bubble[s_i][m_i],
// 			  P[r_s_i][n_m_i],
// 			  {0., beta},
// 			  AbsolutePrecision/U/U,
// 			  true); //it is a bosonic convolution

	  }
	}
      }
    }

    
  
    for( auto r: Bravais_sites ){
      int r_i = create_index(r, L_array);
      for( auto n: Cell_sites ){
	int n_i = create_index(n, std::array{cell_size, cell_size});


	auto dyson_right_hand_side = RPA_Bubble[r_i][n_i] +
	  Int_Bubble_P[r_i][n_i]*U;

	
	  Real norm_max = NormMax( P[r_i][n_i] -
				   dyson_right_hand_side)*U*U;
	difference_between_iterations =
	  std::max(difference_between_iterations, norm_max);

      
	P[r_i][n_i] = P[r_i][n_i]*damping_alpha +
	  dyson_right_hand_side*(1-damping_alpha);
      
      }
    }
    
  
    return difference_between_iterations;

    
  }

}//namespace
