namespace ffd::heat_bath_mc::unit_test{

  template<int i1,
	   int i2,
	   bool is_TI,
	   typename T1,
	   typename T2,
	   typename T3,
	   typename T4>

  auto
  invoke_heat_bath_mc(T1 x1, T2 x2, T3 x3, T4 x4){
    if constexpr(is_TI){
	return HeatBathMC<i1, i2, is_TI>(x1, x2, x3);
      }else{
      return HeatBathMC<i1, i2, is_TI>(x1, x2, x3, x4);
    }
  }
  

}//namespace
