

namespace ffd::user_space{

  template<auto k,
	   template<typename, typename...> typename SetType,
	   typename ValueType,
	   typename... Args>

  auto

  Combination(SetType<ValueType, Args...> set_){
    return ffd::combination::
      Combination<k,
		  SetType,
		  ValueType,
		  Args...>(set_);
  }
  

}//namespace
