namespace ffd::chebyshev_polynomial{

  std::vector<Real>
  ReturnsChebyshevNodesOfOrder1From2To3(int polynomial_order,
					Real lower_limit = -1.,
					Real upper_limit =  1.){
    std::vector<Real> nodes(polynomial_order);
    auto const&a = lower_limit;
    auto const&b = upper_limit;
    for(int j=0; j < polynomial_order; ++j){
      nodes[j] = .5*(a+b) + .5*(b-a)*cos(ffd::core_math::Pi*(j+0.5)/polynomial_order);
    }
    return nodes;
  }

}//namespace
