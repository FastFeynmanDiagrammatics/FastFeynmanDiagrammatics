namespace ffd::user_space{

  template<auto d,
	   template<typename, typename...> typename SetType,
	   typename ValueType,
	   typename... Args>

  auto
  
  CartesianPower(SetType<ValueType, Args...> Set_){
    std::array<SetType<ValueType, Args...>, d> Sets_;
    Sets_.fill(Set_);
    return ffd::cartesian_product::CartesianProduct<d,
						    SetType,
						    ValueType,
						    Args...>(Sets_);
  }
      


}//namespace
