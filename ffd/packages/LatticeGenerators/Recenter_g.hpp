namespace ffd::lattice_g{

  template<int d, int n_atoms = 1>
  auto
  Recenter_g(std::array<int, d> const& L){
    auto ret =
      [L]
      (coord_d<d, n_atoms> const& x){
	auto y = x;
	for(uint j=0; j<d; ++j){
	  while(y[j] < 0) y[j] += L[j];
	  while(y[j] >= L[j]) y[j] -= L[j];
	}
	if constexpr(n_atoms != 1){
	    while(y[d] < 0) y[d] += n_atoms;
	    while(y[d] >= n_atoms) y[d] -= n_atoms;
	  }
	return y;
      };
    return ret;
  }

}//namespace
