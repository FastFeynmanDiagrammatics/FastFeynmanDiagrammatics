namespace ffd::normal_propagator_g{
  
  template<int d, int n_atoms = 1,
  	   class field_d = Real>
  auto
  matrix_diag(std::array<int, d> const& L,
	      ffd::user_space::
	      QFSum<ffd::lattice_g::
	      coord_d<d, n_atoms>, field_d> const& H0){
    using namespace ffd::user_space;
    using sp_d = ffd::lattice_g::coord_d<d, n_atoms>;
    //using st_d = std::pair<Real, sp_d>;
    using ffd::core_math::Pi;
    
    
    //building H_{k,\sigma}
    auto b_oracle =
      BlockOracle_g(H0);
    auto blocks = b_oracle();
    auto ret_index = ffd::index_lattice::
      g<d, n_atoms>(L, b_oracle);
    auto index_f = ret_index.first;
    
    
    //std::size_t counter = 0;
    std::vector<std::vector<std::pair<std::vector<Real>,
				      std::vector<std::vector<Complex>>>>> Eig_k_blocks;
    std::vector<std::array<int, d>> k_array;
    for(uint u=0; u<size(blocks.split)-1; ++u){
      std::vector<std::pair<std::vector<Real>,
			    std::vector<std::vector<Complex>>>>	Eig_k;
      std::size_t const size_block = blocks.split[u+1] - blocks.split[u];
      std::size_t const n_comp = size_block*n_atoms;
      k_array.clear();
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

	    // std::cerr << pos0 << ' ' << pos1 << ':' << k_dx << '\n';
	    H_k[pos0][pos1] += std::exp(Complex(0., k_dx))*value_term;
	  }
	}
	Eig_k.push_back(ffd::eigenvalues_vectors::
			EigenvaluesVectors(H_k));
	k_array.push_back(k);
      }
      Eig_k_blocks.push_back(Eig_k);
    }
    return std::make_tuple(Eig_k_blocks, k_array, b_oracle);
  }

  

}//namespace


