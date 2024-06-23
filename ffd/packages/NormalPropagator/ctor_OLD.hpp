namespace ffd::bravais_propagator{

  template<typename action_t,
	   typename time_t,
	   typename lattice_t>
  auto
  BravaisPropagator(action_t const& S0,
		    time_t Time,
		    lattice_t Lattice,
		    Real precision_G0){
    if(precision_G0 < 0 ){
      precision_G0 = precision_G0_default;
    }
    
    using coord_t = typename
      std::decay_t<decltype(ffd::user_space::CreateCoordinates(Time, Lattice))>;
    using ffd::core_math::Pi;
    //L ??
    L = ffd::lattice::GetLinearSizes(Lattice);
    //get Beta
    beta_ = Beta(Time);
    constexpr int d = ffd::lattice::dimension(Lattice);

    
    //building H_{k,\sigma}
    auto const field_comp = action_fields_components(S0);
    std::size_t const n_comp = size(field_comp);
    std::vector<std::pair<std::vector<Real>,
			  std::vector<std::vector<Complex>>>>
		Eig_k;
    for(auto const& k: ffd::user_space::VectorRange<d>(L)){
      std::vector<std::vector<Complex>>
	H_k(n_comp, std::vector<Complex>(n_comp, 0.));
      
      
      for(uint j=0; j<size(S0.fields); ++j){
	assert((size(S0.fields[j])==2));
	Complex value_term = S0.coef[j];
	std::array<coord_t, 2> positions;
	std::array<ffd::user_space::QF, 2> fields;
	for(int pos: {0, 1}){
	  fields[pos] = S0.fields[j][pos].first;
	  positions[pos] = S0.fields[j][pos].second;
	}

		
	if(fields[0].direction ==
	   ffd::phys::ou){
	  if(fields[0].statistics ==
	     ffd::phys::fermi){
	    value_term *= -1.;
	  }
	  std::swap( fields[0], fields[1] );
	  std::swap( positions[0], positions[1] );
	}
	
	
	assert((fields[0].direction == phys::in));
	assert((fields[1].direction == phys::ou));
	auto dx0 = RealSpace(positions[0]);
	auto dx1 = RealSpace(positions[1]);
	Real k_dx = 0.;
	for(int j=0; j<d; ++j){
	  k_dx += (k[j]*(dx1[j]-dx0[j])*2.)*Pi/L[j];
	}
	auto pos0 = position_in_vector(fields[0], field_comp);
	auto pos1 = position_in_vector(fields[1], field_comp);
	H_k[pos0][pos1] += std::exp(Complex(0., k_dx))*value_term;
      }
      Eig_k.push_back(ffd::eigenvalues_vectors::
		      EigenvaluesVectors(H_k));
    }
    //building the chebyshev polynomials
    for(auto const& x: ffd::user_space::VectorRange<d>(L)){
      for(int s1=0; s1<n_comp; ++s1){
	for(int s2=0; s2<n_comp; ++s2){

	  
	  auto G0_temp =
	    [=](Real tau)->Real{
	      Complex ret = 0.;
	      for(auto k: ffd::user_space::VectorRange<d>(L)){
		auto [index, discard] = compute_index_array<d>(k, L);
		auto e_k = Eig_k[index];
		auto E_k = e_k.first;
		auto Eig_k = e_k.second;
		Complex ret_k = 0.;
		for(uint j=0; j<size(E_k); ++j){
		  ret_k += .5*Eig_k[j][s1]*std::conj(Eig_k[j][s2])*
		    exp(-tau*E_k[j])/cosh(.5*beta_*E_k[j]);
		}
		Real k_x = 0;

		for(int j=0; j<d; ++j){
		  k_x += (k[j]*x[j]*2.)*Pi/L[j];
		}

	      
		ret += std::exp(Complex(0., k_x))*ret_k;
	      }

	      Real ll = 1;
	      for(int j=0; j<d; ++j){
		ll *= L[j];
	      }
	      return std::real(ret/ll);
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
	auto pos1 = position_in_vector(pair.fields[1].first, field_comp);
	auto x1 = pair.fields[0].second;
	auto x2 = pair.fields[1].second;
	auto index = index_difference_coordinates<d>(x1, x2, L);
	auto dt = ffd::get<0>(x1)()-ffd::get<0>(x2)();
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


	return G0[pos0+n_comp*(pos1+n_comp*index)];
      };
    return G0_l;
  }

}//namespace
