namespace ffd::lattice_g{

  template<int d, int n_atoms = 1,
	   typename realspace_f = int>
  auto
  Distance2_g(realspace_f const& Realspace){
    auto Distance =
      [Realspace]
      (coord_d<d, n_atoms> const& x0,
       coord_d<d, n_atoms> const& x1){
	auto const d0 = Realspace(x0);
	auto const d1 = Realspace(x1);
	Real ret = 0.;
	for(uint j=0; j<d; ++j){
	  Real const d01 = std::abs(d0[j]-d1[j]);
	  ret += d01*d01;
	}
	return std::sqrt(ret);
      };
    return Distance;
  }
  
}//namespace
