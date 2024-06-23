namespace ffd::chebyshev_polynomial::unit_test{

  struct easy_to_fit{
    Real a;

    easy_to_fit(Real a_): a(a_) {}
    
    Real operator()(Real x) const{
      Real y = a*x;
      return (1+y+y*y)/cosh(y);
    }
    
  };

}//namespace 
