namespace ffd::lattice_g{

  template<int d, int n_atoms = 1>

  auto
  Realspace_diff_g(std::array<int, d> const& L,
		   std::array<std::array<Real, d>, d> const& bravais,
		   std::array<std::array<Real, d>, n_atoms> const& shift_atoms
		   = std::array<std::array<Real, d>, n_atoms>()){
    auto diff =
      [bravais, shift_atoms]
      (coord_d<d, n_atoms> const& x0,
       coord_d<d, n_atoms> const& x1)
      {
	std::array<Real, d> pos;
	pos.fill(0.);
	for(uint j=0; j<d; ++j){
	  for(uint u=0; u<d; ++u){
	    pos[u] += (x0[j]-x1[j])*
	      bravais[j][u];
	  }
	}
	if constexpr(n_atoms > 1){
	    for(uint u=0; u<d; ++u){
	      pos[u] += shift_atoms[x0[d]][u]
		-shift_atoms[x1[d]][u];
	    }
	  }
	return pos;
      };
    
    
    auto ret =
      [diff, L]
      (coord_d<d, n_atoms> const& x0,
       coord_d<d, n_atoms> const& x1){
	std::array<std::array<int, 2>, d> s_array;
	for(uint j=0; j<d; ++j){
	  s_array[j][0] = -2;
	  s_array[j][1] = 3;
	}
	coord_d<d, n_atoms> x_0, x_1;
	x_1.fill(0);
	if constexpr(n_atoms > 1){
	    x_1[d] = x1[d];
	    while(x_1[d] < 0) x_1[d] += n_atoms;
	    while(x_1[d] >= n_atoms) x_1[d] -= n_atoms;
	    x_0[d] = x0[d];
	    while(x_0[d] < 0) x_0[d] += n_atoms;
	    while(x_0[d] >= n_atoms) x_0[d] -= n_atoms;
	  }
	for(uint j=0; j<d; ++j){
	  int const lbound = -(L[j]-1)/2;
	  auto recenter = ffd::periodic_numbers::
	    Periodic_g<int>({lbound, lbound+L[j]});
	  auto [y, discard] = recenter(x0[j]-x1[j]);
	  x_0[j] = y;
	}
	std::array<Real, d> pos_min;
	Real norm_pos_min = 0.;
	bool first_time = true;
	for(auto s: ffd::user_space::
	      CartesianProductRange(s_array)){
	  auto y_0 = x_0;
	  for(uint j=0; j<d; ++j){
	    y_0[j] += s[j]*L[j];
	  }
	  auto pos = diff(y_0, x_1);
	  Real norm_pos = 0.;
	  for(uint j=0; j<d; ++j){
	    norm_pos += pos[j]*pos[j];
	  }
	  if(first_time || norm_pos < norm_pos_min){
	    pos_min = pos;
	    norm_pos_min = norm_pos;
	    first_time = false;
	  }
	}
	return pos_min;
      };
    
    
    return ret;
  }

}//namespace
