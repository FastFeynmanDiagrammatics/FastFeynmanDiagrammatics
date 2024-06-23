namespace ffd::bravais_propagator{

  template<typename action_t>
  std::vector<ffd::user_space::QF>
  action_fields_components(action_t const& S0){
    std::vector<ffd::user_space::QF> qfields;
    for(auto const& term: S0.fields){
      for(auto const& field_st: term){
	if(position_in_vector(field_st.first,
			      qfields) ==
	   size(qfields) &&
	   field_st.first.direction ==
	   phys::in){
	  qfields.push_back(field_st.first);
	}
      }
    }
    return qfields;
  }
  
}//namespace
