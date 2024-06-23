namespace ffd::chebyshev_f{

  template<class field_d>
  field_d
  Evaluate(std::vector<field_d> const& Coef,
	   std::array<Real, 2> const& bounds,
	   Real x){    
    if(Coef.size() == 0){
      return 0.;
    }
    Real const t = (x - .5*(bounds[0]+bounds[1]))*2./(bounds[1]-bounds[0]);
    std::vector<field_d> d_Cheby(size(Coef)+2, 0.);
    for(int j=size(Coef)-1; j >= 1; j--){
      d_Cheby[j] = 2.*t*d_Cheby[j+1] - d_Cheby[j+2] + Coef[j];
    }
    return t*d_Cheby[1] - d_Cheby[2] + 0.5*Coef[0];
  }

  
  template<class field_d>
  field_d
  Evaluate(std::vector<field_d> const& Coef,
	   std::array<std::size_t, 2> const& start_size,
	   std::array<Real, 2> const& bounds,
	   Real x){    
    if(start_size[1] == 0){
      return 0.;
    }
    Real const t = (x - .5*(bounds[0]+bounds[1]))*2./(bounds[1]-bounds[0]);
    field_d dj2 = 0., dj1 = 0., dj = 0.;
    for(int j=start_size[1]-1; j >= 1; j--){
      dj = 2.*t*dj1 - dj2 + Coef[start_size[0] + j];
      dj2 = dj1;
      dj1 = dj;
    }
    return t*dj1 - dj2 + 0.5*Coef[start_size[0]];
  }


}//namespace
