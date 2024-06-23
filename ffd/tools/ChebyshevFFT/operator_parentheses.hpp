

namespace ffd::chebyshev_fft{

  Real
  ChebyshevFFT::
  operator()(Real x) const{
    auto const [half_mean, half_diff] =
      return_half_mean_diff(LowerUpperBound);
    Real t = (x - half_mean)/half_diff;

    
    std::vector<Real> recursion_vec(size(Coef)+2);
    recursion_vec[size(Coef)+1] = 0;
    recursion_vec[size(Coef)] = 0;
    for(int j=size(Coef)-1; j >= 1; --j){
      recursion_vec[j] = 2*t*recursion_vec[j+1] - recursion_vec[j+2] + Coef[j];
    }
    return t*recursion_vec[1] - recursion_vec[2] + 0.5*Coef[0];
  }


}//namespace
