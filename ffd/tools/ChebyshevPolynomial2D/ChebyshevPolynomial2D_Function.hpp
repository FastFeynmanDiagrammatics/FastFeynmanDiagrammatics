

namespace ffd::chebyshev_polynomial_2d{

  template<typename Field>
  ChebyshevPolynomial2D<Field>::
  ChebyshevPolynomial2D(std::function<Field(Real, Real)>
			function_,
			
			int order_,
						      
			std::array<std::array<Real, 2>, 2>
			rectangular_domain_):
    Order(order_), Domain(rectangular_domain_){
    using ffd::core_math::Pi;

    std::array<std::vector<Real>, 2> nodes;
    for(auto j: {0, 1}){
      auto [lower_bound, upper_bound] = Domain[j];
      nodes[j] = ffd::chebyshev_polynomial::
	ReturnsChebyshevNodesOfOrder1From2To3(Order,
					      lower_bound,
					      upper_bound);
    }


    Coef = std::vector<Field>(Order*Order, 0.);
    std::vector<Real> Cos_k_x(Order*Order);
    for(int x=0; x < Order; ++x){
      for(int k=0; k < Order; ++k){
	Cos_k_x[k+x*Order] = std::cos(Pi*k*(x+.5)/Order);
      }
    }
    for(int r = 0; r < Order*Order; ++r){
      int y = r/Order;
      int x = r%Order;
      Field function_value = function_(nodes[0][x], nodes[1][y]);
      for(int k = 0; k < Order*Order; ++k){
	int ky = k/Order;
	int kx = k%Order;
	Coef[k] += function_value*Cos_k_x[kx+x*Order]*Cos_k_x[ky+y*Order];
      }
    }
    
    for(int k = 0; k < Order*Order; ++k){
      Real factor_normalization = 2./Order;
      Coef[k] *= std::pow(factor_normalization, 2);
    }

  }


}//namespace
