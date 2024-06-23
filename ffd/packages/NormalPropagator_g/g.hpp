namespace ffd::normal_propagator_g{
  
  template<int d, int n_atoms=1,
	   template<class> class chebyshev_t = ffd::chebyshev_polynomial::ChebyshevPolynomial,
  	   class field_d = Real>
  auto
  g(Real Beta,
    std::array<int, d> const& L,
    ffd::user_space::
    QFSum<ffd::lattice_g::
    coord_d<d, n_atoms>, field_d> const& H0,
    Real precision_G0 = 1e-10){
    
    
    using namespace ffd::user_space;
    using sp_d = ffd::lattice_g::coord_d<d, n_atoms>;
    using st_d = std::pair<Real, sp_d>;
    using ffd::core_math::Pi;
    Real const Beta_2 = .5*Beta;
    
    
    //building H_{k,\sigma}
    auto b_oracle =
      BlockOracle_g(H0);
    auto blocks = b_oracle();
    auto ret_index = ffd::index_lattice::
      g<d, n_atoms>(L, b_oracle);
    auto index_f = ret_index.first;
    
    
    VectorContainer<field_d> G0_cheby;
    std::size_t counter = 0;
    for(uint u=0; u<size(blocks.split)-1; ++u){
      std::vector<std::pair<std::vector<Real>,
			    std::vector<std::vector<Complex>>>>	Eig_k;
      std::size_t const size_block = blocks.split[u+1] - blocks.split[u];
      std::size_t const n_comp = size_block*n_atoms;
      std::vector<std::array<int, d>> k_array;
      for(auto const& k: ffd::user_space::VectorRange<d>(L)){
	std::vector<std::vector<Complex>>
	  H_k(n_comp, std::vector<Complex>(n_comp, 0.));
	
	
	for(uint term=0; term<size(H0.fields)/2; ++term){
	  auto [block, pos] = b_oracle(H0.fields[2*term].first);
	  if(block == u){
	    Complex value_term = H0.coef[term];
	    // std::cerr << u << ' '
	    // 	      << term << ' '
	    // 	      << value_term << '\n';
	    std::array<sp_d, 2> x;
	    std::array<QField, 2> qfields;
	    for(int pos: {0, 1}){
	      qfields[pos] = H0.fields[2*term + pos].first;
	      x[pos] = H0.fields[2*term + pos].second;
	    }
	    
	    
	    if(Direction(qfields[1]) ==
	       phys::ou){
	      if(Statistics(qfields[1]) ==
		 phys::fermi){
		value_term *= -1.;
	      }
	      std::swap( qfields[0], qfields[1] );
	      std::swap( x[0], x[1] );
	    }
	    
	    assert((Direction(qfields[1]) == phys::in));
	    assert((Direction(qfields[0]) == phys::ou));
	    Real k_dx = 0.;
	    for(int j=0; j<d; ++j){
	      k_dx += (k[j]*(-x[0][j]+x[1][j])*2.)*Pi/L[j];
	    }

	    auto [discard0, pos0] = b_oracle(qfields[0]);
	    auto [discard1, pos1] = b_oracle(qfields[1]);
	    assert((discard0 == discard1));
	    if constexpr(n_atoms > 1){
		pos0 += size_block*x[0][d];
		pos1 += size_block*x[1][d];
	      }
	    // std::cerr << pos0
	    // 	      << ' '
	    // 	      << ToString(qfields[0])
	    // 	      << x[0]
	    // 	      << ' '
	    // 	      << pos1
	    // 	      << ' '
	    // 	      << ToString(qfields[1])
	    // 	      << x[1]
	    // 	      << ' '
	    // 	      << value_term
	    // 	      << '\n';

	    // std::cerr << pos0 << ' ' << pos1 << ':' << k_dx << '\n';
	    H_k[pos0][pos1] += std::exp(Complex(0., k_dx))*value_term;
	  }
	}
	// for(uint j=0; j<n_comp; ++j){
	// 	for(uint u=0; u<n_comp; ++u){
	// 	  std::cerr << H_k[j][u] << ' ';
	// 	}
	// 	std::cerr << '\n';
	// }
	Eig_k.push_back(ffd::eigenvalues_vectors::
			EigenvaluesVectors(H_k));
	k_array.push_back(k);
      }
      // std::cerr << "eig = " << Eig_k[0].first[0] << ' '
      // 	      << Eig_k[0].first[1] << std::endl;
      
      //building the chebyshev polynomials
      for(std::size_t p1=0; p1<size_block; ++p1){
	for(std::size_t a1=0; a1<n_atoms; ++a1){
	  auto const s1 = p1+size_block*a1;
	  for(std::size_t p0=0; p0<size_block; ++p0){
	    for(std::size_t a0=0; a0<n_atoms; ++a0){
	      auto const s0 = p0+size_block*a0;
	      for(auto const& x0: ffd::user_space::VectorRange<d>(L)){
		auto G0_temp =
		  [Eig_k, k_array, s0, s1, x0, L, Beta, Beta_2](Real tau)
		  -> field_d
		  {
		    Complex ret = 0.;
		    for(std::size_t k_ind=0; k_ind<size(Eig_k); ++k_ind){
		      auto k = k_array[k_ind];
		      auto e_k = Eig_k[k_ind];
		      auto E_k = e_k.first;
		      auto v_k = e_k.second;
		      Complex ret_k = 0.;
		      for(uint j=0; j<size(E_k); ++j){
			//////ONLY WORKS FOR FERMIONS right now
			ret_k += .5*v_k[j][s0]*std::conj(v_k[j][s1])*
			  exp((Beta_2-tau)*E_k[j])/cosh(Beta_2*E_k[j]);
		      }
		      
		      Real k_x = 0;
		      for(int j=0; j<d; ++j){
			k_x += (k[j]*x0[j]*2.)*Pi/L[j];
		      }
		      
		      
		      ret += std::exp(Complex(0., k_x))*ret_k;
		    }
		    
		    std::size_t L_prod = 1;
		    for(int j=0; j<d; ++j) L_prod *= L[j];
		    
		    if constexpr(std::is_same_v<field_d, Real>){
		return std::real(ret/Real(L_prod));
		}else{
		      return ret/Real(L_prod);
		    }
		  };
	    
	    
		auto cheby_temp = 
		  chebyshev_t<field_d>(G0_temp, {0., Beta}, precision_G0);
		cheby_temp.Purify(precision_G0);
		G0_cheby.push_back(cheby_temp.Coef);
		// sp_d x_0, x_1;
		// x_1.fill(0);
		// for(uint j=0; j<d; ++j) x_0[j] = x0[j];
		// if constexpr(n_atoms > 1){
		//     x_0[d] = a0;
		//     x_1[d] = a1;
		//   }
		// std::cerr << "counter = "
		// 	  << counter
		// 	  << ' '
		// 	  << index_f(std::make_pair(blocks[u][p0],
		// 				    x_0),
		// 		     std::make_pair(blocks[u][p1],
		// 				    x_1))
		// 	  << ' '
		// 	  << p0
		// 	  << ' '
		// 	  << x_0
		// 	  << ' '
		// 	  << x_1
		// 	  << '\n';
		counter++;
	      }
	    }
	  }
	}
      }
    }
    
    auto recenter_b = ffd::periodic_numbers::
      Periodic_g({0., Beta}, 1.);
    auto recenter_f = ffd::periodic_numbers::
      Periodic_g({0., Beta}, -1.);

    //    auto evaluate = ffd::chebyshev_f::Evaluate_g<field_d>({0., Beta});

    auto G0 =
      [G0_cheby, index_f, Beta, recenter_f, recenter_b
       //,evaluate
       ]
      (std::pair<QField, st_d> const& QX0,
       std::pair<QField, st_d> const& QX1)
      -> field_d
      {
	auto [q0, X0] = QX0;
	auto [q1, X1] = QX1;
	assert((Direction(q0) !=
		Direction(q1)));
	Real sign_swap = 1.;
	if(Direction(q0) == phys::ou){
	  if(Statistics(q0) == phys::fermi) sign_swap = -1.;
	  std::swap(q0, q1);
	  std::swap(X0, X1);
	}

	// Real phase_itime = 1.;
	// assert((Statistics(q0) ==
	// 	Statistics(q1)));
	// if(Statistics(q0) == phys::fermi) phase_itime = -1.;
	
	
	auto const [t0, x0] = X0;
	auto const [t1, x1] = X1;
	Real const eps = 10*std::numeric_limits<Real>::epsilon()*
	  1;
	  // (1-2*(IsNambu(q0) && ffd::user_space::q_field::
	  // 	value_bit(q0,
	  // 		  user_space::q_field::sign_bit)));
	auto const [dt, sign_time] = ((Statistics(q0) == phys::fermi) ?
				      recenter_f(t0-t1-eps) :
				      recenter_b(t0-t1-eps));

	//////////TESTTTTTTTTTTTTTTTT////////////
	// sign *= index_f(std::make_pair(FlipSpin(q0), x1),
	// 		   std::make_pair(FlipSpin(q1), x0));
	// sign *= index_f(std::make_pair(q0, x1),
	// 		std::make_pair(q1, x0));
	///////////////TESTTTTTTTTTTTTTT
      
	return sign_time*sign_swap*
	  ffd::chebyshev_f::Evaluate(G0_cheby.container,
				     G0_cheby.start_size(index_f(std::make_pair(q0, x0),
								 std::make_pair(q1, x1))),
				     {0., Beta},
				     dt);
      };
    
    
    return std::make_pair(G0, b_oracle);
  }
  
}//namespace


