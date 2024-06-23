

namespace ffd::chebyshev_polynomial::unit_test{

  template<typename T1, typename T2>
  Real compute_max_error(T1 const& func1,
			 T2 const& func2,
			 Real a = -1, Real b = 1){
    const long N_points = 1<<16;
    Real ret = 0.;
    for(int j=0; j <= N_points; ++j){
      Real x = a + (b-a)*j*1./N_points;
      Real diff = std::abs(func1(x)-func2(x));
      ret = diff > ret ? diff : ret;
    }
    return ret;
  }

}//namespace
