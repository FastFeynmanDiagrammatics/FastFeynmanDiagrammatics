namespace ffd::imaginary_time_convolution_2d{


  template<typename f_t,
	   typename g_t,
	   typename h_t>

  
  dflt::chebyshev_2d_t<Real>
  
  VeryFastFermionicConvolute123FromZeroTo4Dflt(f_t&& f_v,
					       g_t&& g_v,
					       h_t&& h_v,
					       Real beta,
					       Real precision){
    return VeryFastFermionicConvolute123FromZeroTo4<dflt::chebyshev_2d_t<Real>>(std::forward<f_t>(f_v),
										std::forward<g_t>(g_v),
										std::forward<h_t>(h_v),
										beta,
										precision);
  }




}//namespace
