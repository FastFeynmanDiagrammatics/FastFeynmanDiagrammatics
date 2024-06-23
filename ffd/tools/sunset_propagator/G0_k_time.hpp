namespace ffd::sunset_propagator{

  template<std::size_t l_x,
	   std::size_t l_y,
	   std::size_t n_t,
	   template<typename, std::size_t> typename array_t = ffd::l_array>
  
  array_t<fermi_cheby_s<n_t>, (1ul<<(l_x+l_y))>

  G0_k_time(Real const beta,
	    Real const mu0,
	    Real const t = 1.,
	    Real const tp = 0.){
    using ffd::core_math::Pi;
    using std::cos, std::exp, std::cosh;
    

    std::size_t constexpr size_x_y = (1ul<<(l_x+l_y));
    std::size_t constexpr  size_x =  (1ul<<l_x);
    std::size_t constexpr  size_y =  (1ul<<l_y);

    
    array_t<fermi_cheby_s<n_t>, size_x_y> G0_kt;
    for(std::size_t ny=0; ny<size_y; ++ny){
      for(std::size_t nx=0; nx<size_x; ++nx){
	Real const kx = Real(nx)*2.*Pi/size_x;
	Real const ky = Real(ny)*2.*Pi/size_y;
	Real const cos_kx = cos(kx);
	Real const cos_ky = cos(ky);
	Real const xi_k = -2.*t*(cos_kx+cos_ky) - 4.*tp*cos_kx*cos_ky - mu0;

	
	G0_kt[nx+size_x*ny] =
	  fermi_cheby_s<n_t>([=](Real const tau){
			       return -.5*exp((.5*beta-tau)*xi_k)
				 /cosh(.5*beta*xi_k);}, {0., beta});
      }
    }
    
    
    return G0_kt;
  }
    
}//namespace
