

namespace ffd::chebyshev_polynomial::unit_test{

  struct not_so_easy_to_fit{
    Real a, b;

    not_so_easy_to_fit(int a_, int b_): a(a_), b(b_) {}
        
    Real operator()(Real x) const{
      return std::sin(a*sin(b*x));
    }
    
  };
  
}
