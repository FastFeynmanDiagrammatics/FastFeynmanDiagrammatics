namespace ffd::gauss_kronrod{

  template<typename func_t>
  Real
  Integrate(func_t&& func_v,
	    std::array<Real, 2> interval,
	    Real precision){
    auto [value3, error3] = Integrate_s<3>(func_v, interval);
    if( error3 < precision ){
      return value3;
    }


    // auto [value7, error7] = Integrate_s<7>(func_v, interval);
    // if( error7 < precision ){
    //   return value7;
    // }
    
    
    
    auto [value15, error15] = Integrate_s<15>(func_v, interval);
    if( error15 < precision ){
      return value15;
    }



    // auto [value31, error31] = Integrate_s<31>(func_v, interval);
    // if( error31 < precision ){
    //   return value31;
    // }

    
    auto [value63, error63] = Integrate_s<63>(std::forward<func_t>(func_v), interval);

    
#ifndef FFD_GAUSSKRONROD_DO_NOT_CHECK_ERROR
    assert(( error63 < precision ));
#endif

    
    return value63;
  }
  
}//namespace
