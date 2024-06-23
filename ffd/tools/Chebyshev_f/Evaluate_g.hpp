namespace ffd::chebyshev_f{

  template<class field_d>
  auto
  Evaluate_g(std::array<Real, 2> const& bounds){
    Real const mean = .5*(bounds[0]+bounds[1]);
    Real const norm = 2./(bounds[1]-bounds[0]);

    auto ret =
      [mean, norm]
      (std::vector<field_d> const& Coef,
       std::array<std::size_t, 2> const& start_size,
       Real x)
      {
	if(start_size[1] == 0){
	  return 0.;
	}
	Real const t = (x - mean)*norm;
	field_d dj2 = 0., dj1 = 0., dj = 0.;
	for(int j=start_size[1]-1; j >= 1; j--){
	  dj = 2.*t*dj1 - dj2 + Coef[start_size[0] + j];
	  dj2 = dj1;
	  dj1 = dj;
	}
	return t*dj1 - dj2 + 0.5*Coef[start_size[0]];
      };

    return ret;
  }


}//namespace
