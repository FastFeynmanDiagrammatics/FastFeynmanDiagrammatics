

namespace ffd::chebyshev_polynomial::unit_test{
  
  std::vector<Real> compute_vector_of_values(std::function<Real(Real)> func_, int N_Cheby, Real a = -1., Real b = 1.){
    std::vector<Real> ret;
    for(int j=0; j < N_Cheby; ++j){
      ret[j] = func_(.5*(a+b)+.5*(b-a)*cos(ffd::core_math::Pi*(j+.5)/N_Cheby));
    }
    return ret;
  }

}//namespace
