namespace ffd::user_space{

  template<auto... sequence,
	   typename value_t>

  auto
  CreateSafeMap(std::array<value_t, sizeof...(sequence)> values_){
    return SafeMap(IntSequence<sequence...>(), values_);
  }


}//namespace
