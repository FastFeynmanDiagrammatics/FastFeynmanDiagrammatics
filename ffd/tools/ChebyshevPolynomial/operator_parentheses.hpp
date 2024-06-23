namespace ffd::chebyshev_polynomial{

  template<typename Field>
  Field ChebyshevPolynomial<Field>::operator()(Real x) const{
    if(Coef.size() == 0){
      return 0.;
    }
    Real t = (x - .5*(LowerLimit+UpperLimit))*2./(UpperLimit-LowerLimit);
    Field dj2 = 0., dj1 = 0., dj = 0.;
    for(int j=Coef.size()-1; j >= 1; j--){
      dj = 2.*t*dj1 - dj2 + Coef[j];
      dj2 = dj1;
      dj1 = dj;
    }
    return t*dj1 - dj2 + 0.5*Coef[0];
  }

}//namespace
