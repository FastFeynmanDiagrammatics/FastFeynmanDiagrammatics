namespace ffd::integrate{

  template<Method method_name,
	   typename function_t>
  Real Integrate(function_t F,
		 std::array<Real, 2> interval,
		 Real AbsolutePrecision){
    if constexpr(method_name == Method::ClenshawCurtis){
	return Integrate_clenshaw(F, interval, AbsolutePrecision);
      }else{
      assert(false);
    }
  }


}//namespace
