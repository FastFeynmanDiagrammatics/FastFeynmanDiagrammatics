namespace ffd::user_space{

  template<auto d,
	   template<typename, typename...> typename SetType,
	   typename ValueType,
	   typename... Args>

  auto

  CartesianProduct(std::array<SetType<ValueType, Args...>, d> Sets_){
    return ffd::cartesian_product::CartesianProduct<d,
						    SetType,
						    ValueType,
						    Args...>(Sets_);
  }

  

  template<auto d,
	   template<typename, typename...> typename SetType,
	   typename ValueType,
	   typename... Args>

  auto

  CartesianProduct(SetType<ValueType, Args...> Set_){
    std::array<SetType<ValueType, Args...>, d> Sets;
    Sets.fill(Set_);

    
    return ffd::cartesian_product::
      CartesianProduct<d,
		       SetType,
		       ValueType,
		       Args...>(Sets);
  }


}//namespace
