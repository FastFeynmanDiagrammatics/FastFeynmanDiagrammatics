namespace ffd::lattice_g{

  template<int d, int n_atoms = 1>

  auto
  Realspace_unary_g(std::array<int, d> const& L,
		    std::array<std::array<Real, d>, d> const& bravais,
		    std::array<std::array<Real, d>, n_atoms> const& shift_atoms
		    = std::array<std::array<Real, d>, n_atoms>()){
    auto ret =
      [L, bravais, shift_atoms]
      (coord_d<d, n_atoms> const& x){
	using namespace ffd::user_space;
	std::array<Real, d> pos;
	pos.fill(0.);
	for(uint j=0; j<d; ++j){
	  int const lbound = -(L[j]-1)/2;
	  auto recenter = ffd::periodic_numbers::
	    Periodic_g<int>({lbound, lbound+L[j]});
	  auto [y, discard] = recenter(x[j]);
	  for(uint u=0; u<d; ++u){
	    pos[u] += y*bravais[j][u];
	  }
	}
	if constexpr(n_atoms > 1){
	    auto x_d = x[d];
	    while(x_d<0) x_d += n_atoms;
	    while(x_d>=n_atoms) x_d -= n_atoms;
	    for(uint j=0; j<d; ++j){
	      pos[j] += shift_atoms[x_d][j];
	    }
	  }
	return pos;
      };
    return ret;
  }


}//namespace
