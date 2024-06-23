namespace ffd::gauss_kronrod{

  template<typename func_t>
  Real
  Integrate(func_t&& func_v,
	    std::array<Real, 2> interval,
	    Real AbsolutePrecision = 1e-10);

  
}//namespace
