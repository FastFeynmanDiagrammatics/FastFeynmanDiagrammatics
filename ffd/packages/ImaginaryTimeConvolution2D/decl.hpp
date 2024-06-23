namespace ffd::imaginary_time_convolution_2d{

  namespace dflt{
    template<typename Field>
    using chebyshev_2d_t = ffd::chebyshev_polynomial_2d::ChebyshevPolynomial2D<Field>;
  }
  

  
  template<typename Field,
	   typename chebyshev_2d_t = dflt::chebyshev_2d_t<Field>>
  chebyshev_2d_t
  Convolute123FromZeroTo4(std::function<Field(Real)> f,
			  std::function<Field(Real)> g,
			  std::function<Field(Real)> h,
			  Real beta_,
			  std::array<bool, 2> g_is_fermion__h_is_fermion_,
			  Real absolute_precision_ = 1e-10);



  
  template<typename Field,
	   typename chebyshev_2d_t = 
	   ffd::chebyshev_polynomial_2d::ChebyshevPolynomial2D<Field>>

  chebyshev_2d_t
  
  FastFermionicConvolute123FromZeroTo4(std::function<Field(Real)> f,
				       std::function<Field(Real)> g,
				       std::function<Field(Real)> h,
				       Real beta_,
				       Real absolute_precision_ = 1e-10);


  template<typename chebyshev_2d_t,
	   typename f_t,
	   typename g_t,
	   typename h_t>

  chebyshev_2d_t

  VeryFastFermionicConvolute123FromZeroTo4(f_t&& f_v,
					   g_t&& g_v,
					   h_t&& h_v,
					   Real beta,
					   Real precision);


  template<typename f_t,
	   typename g_t,
	   typename h_t>

  
  dflt::chebyshev_2d_t<Real>
  
  VeryFastFermionicConvolute123FromZeroTo4Dflt(f_t&& f_v,
					       g_t&& g_v,
					       h_t&& h_v,
					       Real beta,
					       Real precision);
  
  
	
}//namespace
