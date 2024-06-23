namespace ffd::user_space{
  
  template<auto... sequence>
  struct IntSequence{
    IntSequence(){};
  };



  namespace safe_map{

    template<auto value, int rec, auto first, auto... sequence>
    constexpr int
    static_position_of_value_in_sequence_rec(){
      if constexpr(value == first){
	  return rec;
	}else if constexpr(sizeof...(sequence) == 0){
	return -1;
      }else{
	return static_position_of_value_in_sequence_rec<value, rec+1, sequence...>();
      }
    }


    template<auto value, auto... sequence>
    constexpr int
    static_position_of_value_in_sequence(){
      return static_position_of_value_in_sequence_rec<value, 0, sequence...>();
    }

    
    template<typename value_t, int rec, auto first, auto... sequence>
    int
    dynamic_position_of_value_in_sequence_rec(value_t value){
      if (value == first){
	  return rec;
      }else{
	if constexpr(sizeof...(sequence) == 0){
	    return -1;
	  }else{
	  return dynamic_position_of_value_in_sequence_rec<value_t, rec+1, sequence...>(value);
	}
      }
    }


    template<typename value_t, auto... sequence>
    int
    dynamic_position_of_value_in_sequence(value_t value){
      return dynamic_position_of_value_in_sequence_rec<value_t, 0, sequence...>(value);
    }
    
  }//namespace

}//namespace
