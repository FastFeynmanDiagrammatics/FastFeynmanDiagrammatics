namespace ffd::periodic_numbers{

  template<typename field1_d = Real,
	   typename field2_d = int>
  auto
  Periodic_g(std::array<field1_d, 2> const& bounds,
	     field2_d phase = 1){
    auto const diff = bounds[1]-bounds[0];
    auto ret =
      [bounds, diff, phase]
      (field1_d x){
	field2_d phase_ret = 1;
	while(x < bounds[0]){
	  x += diff;
	  phase_ret *= phase;
	}
	if constexpr(std::is_integral_v<field1_d>){
	while(x >= bounds[1]){
	  x -= diff;
	  phase_ret *= phase;
	}
	}else{
	  while(x > bounds[1]){
	    x -= diff;
	    phase_ret *= phase;
	  }
	}
	return std::make_pair(x, phase_ret);
      };
    return ret;
  }

}//namespace
