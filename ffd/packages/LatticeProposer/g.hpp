namespace ffd::lattice_proposer{

  template<int d, int n_atoms = 1,
	   typename proposer_f=int,
	   typename distance_f=int>
  auto
  g(std::array<int, d> const& L,
    proposer_f const& Proposer,
    distance_f const& Distance){
 
    using ffd::user_space::Proba;  
    
    std::size_t const L_prod = [L](){uint ret=1; for(uint j=0; j<d; ++j) ret *= L[j];return ret;}();
    std::size_t const size_vec = L_prod*n_atoms;
    std::array<std::vector<Real>, n_atoms> proba;
    std::array<std::vector<Real>, n_atoms> cumul;
    std::vector<coord_d<d, n_atoms>> indexes(size_vec);
    auto lattice_index =
      [L](auto const& x){
	std::size_t index = 0, l_prod = n_atoms;
	if constexpr(n_atoms != 1){ index = x[d];}
	for(uint j=0; j<d; ++j){
	  index += x[j]*l_prod;
	  l_prod *= L[j];
	}
	return index;
      };
    for(uint a0=0; a0<n_atoms; ++a0){
      coord_d<d, n_atoms> x0;
      x0.fill(0);
      if constexpr(n_atoms != 1){
	  x0[d] = a0;
	}
      std::vector<Real> vec(size_vec);
      for(uint a1=0; a1<n_atoms; ++a1){
	coord_d<d, n_atoms> x1;
	x1.fill(0);
	if constexpr(n_atoms != 1){
	  x1[d] = a1;
	}
	for(auto const& x: ffd::user_space::VectorRange<d>(L)){
	  for(uint j=0; j<d; ++j) x1[j] = x[j];
	  std::size_t index = lattice_index(x1);
	  indexes[index] = x1;
	  {
	    using namespace ffd::user_space;
	    // std::cerr << index << " indexes = " << indexes[index] << '\n';
	  }
	  vec[index] = Proposer(Distance(x0, x1));
	}
      }
      Real sum = 0.;
      for(uint j=0; j<size(vec); ++j) sum += vec[j];
      for(uint j=0; j<size(vec); ++j) vec[j] /= sum;
      std::vector<Real> cml(size_vec, 0.);
      cml[0] = vec[0];
      for(uint j=1; j<size(vec); ++j) cml[j] = cml[j-1] + vec[j];
      proba[a0] = vec;
      cumul[a0] = cml;
    }
    auto Recenter = ffd::lattice_g::
      Recenter_g<d, n_atoms>(L);
    
    // for(uint j=0; j<size(indexes); ++j){
    //   using namespace ffd::user_space;
    //   std::cerr << j << ' ' << indexes[j]
    //   		<< proba[0][j] << ' '
    //   		<< proba[1][j] << '\n';
    // }
    
    auto conditional =
      [proba, Recenter, lattice_index]
      (coord_d<d, n_atoms> const& x1,
       coord_d<d, n_atoms> const& x0){
	auto x01 = x0;
	for(uint j=0; j<d; ++j){
	  x01[j] -= x1[j];
	}
	x01 = Recenter(x01);
	if constexpr(n_atoms != 1) x01[d] = x0[d];
	auto const index = lattice_index(x01);
	std::size_t a0 = 0;
	if constexpr(n_atoms != 1) a0 = x1[d];
	return proba[a0][index];
      };
    
    
    auto propose =
      [cumul, indexes, Recenter]
      (coord_d<d, n_atoms> const& x0){
	std::size_t index_0 = 0;
	if constexpr(n_atoms != 1) index_0 = x0[d];
	Real const p = Proba();
	auto it = std::lower_bound(cbegin(cumul[index_0]),
				   cend(cumul[index_0]),
				   p);
	std::size_t index = std::distance(begin(cumul[index_0]),
					  it);
	coord_d<d, n_atoms> x1;
	for(uint j=0; j<d; ++j){
	  x1[j] = x0[j] + indexes[index][j];
	}
	if constexpr(n_atoms!=1) x1[d] = indexes[index][d];
	return Recenter(x1);
      };


    coord_d<d, n_atoms> origin;
    origin.fill(0);

    
    auto ret =
      [conditional, propose, origin]
      (auto&& ...x){
	if constexpr(sizeof...(x) == 0){
	    return origin;
	  }else{
	  if constexpr(sizeof...(x) == 1){
	      return propose(std::forward<decltype(x)>(x)...);
	    }else{
	    return conditional(std::forward<decltype(x)>(x)...);
	  }
	}
      };
    
    return ret;
  }

}//namespace
