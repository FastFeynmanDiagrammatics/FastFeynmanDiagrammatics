namespace ffd::bravais_propagator{

  template<typename chebyshev_t = ffd::chebyshev_polynomial::ChebyshevPolynomial<Real>,
	   typename action_t=int,
	   typename time_t=int,
	   typename lattice_t=int>
  auto
  Propagator(action_t const& S0,
	     time_t Time,
	     lattice_t Lattice,
	     Real precision_G0 = 1e-8){
    using coord_t = typename
      std::decay_t<decltype(ffd::user_space::CreateCoordinates(Time, Lattice))>;
    using ffd::core_math::Pi;
    constexpr int d = ffd::lattice::dimension(Lattice);
    //L ??
    auto L = ffd::lattice::GetLinearSizes(Lattice);
    int L_prod = 1;
    for(uint j=0; j<d; ++j) L_prod *= L[j]; 
    //get Beta
    Real const beta_ = Beta(Time);
    Real const beta_2 = .5*beta_;
    
        
    //building H_{k,\sigma}
    auto const field_comp = action_fields_components(S0);
    
    // for(uint j=0; j<size(field_comp); ++j){
    //   std::cerr << field_comp[j] << ' ';
    // }
    // std::cerr << "---------------\n";

    std::size_t const n_comp = size(field_comp);
    std::vector<std::pair<std::vector<Real>,
			  std::vector<std::vector<Complex>>>>
		Eig_k;
    for(auto const& k: ffd::user_space::VectorRange<d>(L)){
      std::vector<std::vector<Complex>>
	H_k(n_comp, std::vector<Complex>(n_comp, 0.));
      
      
      for(uint j=0; j<size(S0); ++j){
	assert((size(S0.fields[j])==2));
	Complex value_term = S0.coef[j];
	//	std::cerr << j << ' ' << value_term << '\n';
	std::array<coord_t, 2> positions;
	std::array<ffd::user_space::QF, 2> qfields;
	for(int pos: {0, 1}){
	  qfields[pos] = S0.fields[j][pos].first;
	  positions[pos] = S0.fields[j][pos].second;
	  //	  std::cerr << pos << ' ' << fields[pos] << '\n';
	}
	
	
	if(qfields[0].direction ==
	   ffd::phys::ou){
	  if(qfields[0].statistics ==
	     ffd::phys::fermi){
	    value_term *= -1.;
	  }
	  std::swap( qfields[0], qfields[1] );
	  std::swap( positions[0], positions[1] );
	}
	
	
	assert((qfields[0].direction == phys::in));
	assert((qfields[1].direction == phys::ou));
	auto x0 = ffd::imaginary_time_lattice::
	  GetSpaceCoordinates<lattice_t>(positions[0]);
	auto x1 = ffd::imaginary_time_lattice::
	  GetSpaceCoordinates<lattice_t>(positions[1]);
	Real k_dx = 0.;
	for(int j=0; j<d; ++j){
	  k_dx += (k[j]*(x0[j]-x1[j])*2.)*Pi/L[j];
	}
	auto pos0 = position_in_vector(qfields[0], field_comp);
	auto pos1 = position_in_vector(Bar(qfields[1]), field_comp);
	// std::cerr << pos0 << ' ' << pos1 << ':' << k_dx << '\n';
	H_k[pos0][pos1] += std::exp(Complex(0., k_dx))*value_term;
      }
      // for(uint j=0; j<n_comp; ++j){
      // 	for(uint u=0; u<n_comp; ++u){
      // 	  std::cerr << H_k[j][u] << ' ';
      // 	}
      // 	std::cerr << '\n';
      // }
      Eig_k.push_back(ffd::eigenvalues_vectors::
		      EigenvaluesVectors(H_k));
    }
    // std::cerr << "eig = " << Eig_k[0].first[0] << ' '
    // 	      << Eig_k[0].first[1] << std::endl;
    
    std::vector<chebyshev_t> G0(L_prod*n_comp*n_comp);
    //building the chebyshev polynomials
    for(auto const& x: ffd::user_space::VectorRange<d>(L)){
      for(std::size_t s1=0; s1<n_comp; ++s1){
	for(std::size_t s2=0; s2<n_comp; ++s2){
	  	  
	  auto G0_temp =
	    [=](Real tau)->Real{
	      Complex ret = 0.;
	      for(auto k: ffd::user_space::VectorRange<d>(L)){
		auto [index, discard] = compute_index_array<d>(k, L);
		auto e_k = Eig_k[index];
		auto E_k = e_k.first;
		auto v_k = e_k.second;
		Complex ret_k = 0.;
		for(uint j=0; j<size(E_k); ++j){
		  //////ONLY WORKS FOR FERMIONS right now
		  ret_k += .5*v_k[j][s1]*std::conj(v_k[j][s2])*
		    exp((beta_2-tau)*E_k[j])/cosh(beta_2*E_k[j]);
		}
		Real k_x = 0;
		
		for(int j=0; j<d; ++j){
		  k_x += (k[j]*x[j]*2.)*Pi/L[j];
		}

	      
		ret += std::exp(Complex(0., k_x))*ret_k;
	      }


	      return std::real(ret/Real(L_prod));
	    };
	  
	  
	  auto [index, discard] = compute_index_array<d>(x, L);
	  G0[s1+n_comp*(s2+n_comp*index)] = chebyshev_t(G0_temp, {0., beta_}, precision_G0);
	}
      }
    }

    
    auto G0_l =
      [=](ffd::user_space::
	  QuantumFieldPair<coord_t> const& pair){
	auto pos0 = position_in_vector(pair.fields[0].first, field_comp);
	auto pos1 = position_in_vector(Bar(pair.fields[1].first), field_comp);
	if(std::max(pos0, pos1) == size(field_comp)) return 0.;
	auto x1 = pair.fields[0].second;
	auto x2 = pair.fields[1].second;
	auto index = index_difference_coordinates<d, lattice_t>(x1, x2, L);
	auto dt = ffd::get<0>(x1)() - ffd::get<0>(x2)();
	int sign = 1;
	while(dt<0){
	  if(pair.fields[0].first.statistics ==
	     phys::fermi){
	    sign = -sign;
	  }
	  dt+=beta_;
	}
	while(dt>beta_){
	  if(pair.fields[0].first.statistics ==
	     phys::fermi){
	    sign = -sign;
	  }
	  dt -= beta_;
	}


	return pair.sign*sign*G0[pos0+n_comp*(pos1+n_comp*index)](dt);
      };
    return G0_l;
  }

}//namespace
