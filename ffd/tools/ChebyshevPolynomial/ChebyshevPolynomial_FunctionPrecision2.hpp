

namespace ffd::chebyshev_polynomial{
  
  template<typename Field>
  ChebyshevPolynomial<Field>::
  ChebyshevPolynomial(std::function<Field(Real)> func,
		      std::array<Real, 2> a,
		      Real a2,
		      int a3,
		      Real a4,
		      Real a5):
    ChebyshevPolynomial(func, a[0], a[1], a2, a3, a4, a5){}

}//namespace
