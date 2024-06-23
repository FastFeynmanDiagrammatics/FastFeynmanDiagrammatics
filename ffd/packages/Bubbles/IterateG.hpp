namespace ffd::bubbles{


  template<typename ImaginaryTime_t,
	   typename Lattice_t>

  typename Bubbles<ImaginaryTime_t,
		   Lattice_t>::error_t
  
  Bubbles<ImaginaryTime_t,
	  Lattice_t>::
  IterateG(){
    Real difference_between_iterations = 0.;

    
    using ffd::vector_range::Range;
    using std::size, std::abs;


    GF_t G0_Sigma, G0_Sigma_G;
    auto Bravais_sites = ffd::user_space::
      CartesianProductRange(L_array);
    auto Cell_sites = ffd::user_space::
      CartesianProductRange<2>(cell_size);


    std::array cell_sizes{cell_size, cell_size};
    for( int spin: {0, 1} ){
      G0_Sigma[spin].resize(L_product);
      G0_Sigma_G[spin].resize(L_product);
      for( auto r: Bravais_sites ){
	int r_i = create_index(r, L_array);
	
	for( auto n: Cell_sites ){
	  int n_i = create_index(n, cell_sizes);
	  G0_Sigma[spin][r_i][n_i] = chebyshev_t();
	  G0_Sigma_G[spin][r_i][n_i] = chebyshev_t();
	  
	  for( auto s: Bravais_sites ){
	    auto r_s = substraction_indices<d>(r,
					     s,
					     L_array);

	  int s_i = create_index(s, L_array);
	  int r_s_i = create_index(r_s, L_array);


	    for( auto m: Cell_sites ){


	      auto n_m = substraction_indices<2>(n,
						 m,
						 cell_sizes);
	      int m_i = create_index(m, cell_sizes);
	      int n_m_i = create_index(n_m, cell_sizes);

	    
	      auto Sigma = [r_s_i, n_m_i, spin, this](Real tau)->Real{
			     Real U2P = U*U*P[r_s_i][n_m_i](tau);
			     if(RPA == RPA_scheme::ph){
			       return U2P*G[1-spin][r_s_i][n_m_i](tau);
			     }else{
			       return -U2P*G[1-spin][r_s_i][n_m_i](beta-tau);
			     }
			   };

	      
	      G0_Sigma[spin][r_i][n_i] += // ffd::convolution::
		// Convolution(G0[spin][s_i][m_i],
		// 	    Sigma,
		// 	    {0., beta},				
		// 	    AbsolutePrecision,
		// 	    false); //it is NOT a bosonic convolution
		ffd::imaginary_time_convolution::
		Convolute1And2FromZeroTo3<Real,
					  chebyshev_t>(G0[spin][s_i][m_i],
						       Sigma,
						       beta,
						       true, //fermionic
						       AbsolutePrecision);
	    }
	  }
	}
      }
    }


    

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
	    for( int spin: {0, 1} ){
	      std::array cell_sizes{cell_size, cell_size};
	      auto n_m = substraction_indices<2>(n,
						 m,
						 cell_sizes);
	      int n_i = create_index(n, cell_sizes);
	      int m_i = create_index(m, cell_sizes);
	      int n_m_i = create_index(n_m, cell_sizes);
	    
	  
	      G0_Sigma_G[spin][r_i][n_i] += ffd::imaginary_time_convolution::
		Convolute1And2FromZeroTo3<Real,
					  chebyshev_t>(G0_Sigma[spin][s_i][m_i],
						       G[spin][r_s_i][n_m_i],
						       beta,
						       true, //fermionic
						       AbsolutePrecision);
		// ffd::convolution::
		// Convolution(G0_Sigma[spin][s_i][m_i],
		// 	    G[spin][r_s_i][n_m_i],
		// 	    {0., beta},
		// 	    AbsolutePrecision,
		// 	    false); //it is NOT a bosonic convolution
			    
	      
	    }
	  }
	}
      }
    }


    
    for( int spin: {0, 1} ){    
      for( auto r: Bravais_sites ){
	int r_i = create_index(r, L_array);
	for( auto n: Cell_sites ){
	  int n_i = create_index(n, std::array{cell_size, cell_size});


	  auto dyson_right_hand_side = G0[spin][r_i][n_i] +
	    G0_Sigma_G[spin][r_i][n_i];

	  
	  Real norm_max = NormMax( G[spin][r_i][n_i] -
				   dyson_right_hand_side);
	  difference_between_iterations =
	    std::max(difference_between_iterations, norm_max);

      
	  G[spin][r_i][n_i] = G[spin][r_i][n_i]*damping_alpha +
	    dyson_right_hand_side*(1-damping_alpha);
	}
      }
    }
    
  
    return difference_between_iterations;
    
  }

}//namespace
