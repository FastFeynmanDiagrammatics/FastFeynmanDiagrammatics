

namespace ffd::chebyshev_fft{

  auto
  return_half_mean_diff(std::array<Real, 2> x){
    std::array<Real, 2> mean_diff;
    mean_diff.fill(0.);

    
    for(int j: {0, 1}){
      for(int u: {0, 1}){
	mean_diff[j] += .5*(1-2*j*u)*x[1-u];
      }
    }

    
    return mean_diff;
  }

}//namespace
