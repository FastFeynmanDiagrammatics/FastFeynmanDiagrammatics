namespace ffd::user_space{

  
  template<int component_number,
	   int d,
	   typename... args>
  auto&
  component(ffd::periodic_coordinate::PeriodicCoordinates<d, args...>& C_){
    return ffd::get<component_number>(C_).Variable;
  }


  
  template<int component_number,
	   int d,
	   typename... args>
  auto const&
  component(ffd::periodic_coordinate::PeriodicCoordinates<d, args...> const& C_){
    return ffd::get<component_number>(C_).Variable;
  }


}//namespace
