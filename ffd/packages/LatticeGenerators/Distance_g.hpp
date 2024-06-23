namespace ffd::lattice_g{

  template<int d, int n_atoms = 1,
	   typename realspace_f = int>
  auto
  Distance_g(realspace_f const& Realspace){
    auto Distance =
      [Realspace]
      (coord_d<d, n_atoms> const& x0,
       coord_d<d, n_atoms> const& x1){
	auto const diff = Realspace(x0, x1);
	Real ret = 0.;
	for(uint j=0; j<d; ++j){
	  Real const d01 = std::abs(diff[j]);
	  ret += d01*d01;
	}
	return std::sqrt(ret);
      };
    return Distance;
  }
  
}//namespace
