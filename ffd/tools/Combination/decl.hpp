

namespace ffd::combination{

  template<auto k,
	   template<typename, typename...> typename SetType,
	   typename ValueType,
	   typename... Args>
  
  struct Combination{
    SetType<ValueType, Args...> Set;

    std::array<ValueType, k> Values;
    
    std::array<std::size_t, k> indexes;
    
    

    Combination(SetType<ValueType, Args...> set_){
      for( int j: ffd::vector_range::Range(k) ){
	Set = set_;
	indexes[j] = 0;
      }
    }
    
    
    Combination(SetType<ValueType, Args...> set_,
		std::array<std::size_t, k> indexes_):
      Combination(set_){ indexes = indexes_; }
    

    

    
    Combination operator++(){
      for( int j: ffd::vector_range::Range(k) ){
	if( j == k-1){
	  ++indexes[j];
	}else{
	  if( indexes[j] < indexes[j+1] ){
	    ++indexes[j];
	    break;
	  }else{
	    indexes[j] = 0;
	  }
	}
      }
      return *this;
    }


    bool operator!=(Combination const& it) const{
      return indexes != it.indexes;
    }


    std::array<ValueType, k>& operator*(){
      for( int j: ffd::vector_range::Range(k) ){
	auto value_iterator = Set.begin();
	for( [[maybe_unused]] int iter:
	       ffd::vector_range::Range(indexes[j]) ){
	  ++value_iterator;
	}
	Values[j] = *value_iterator;
      }
      return Values;
    }

    
    Combination begin() const{
      return Combination(Set);
    }


    Combination end() const{
      std::array<std::size_t, k> indexes_end;
      indexes_end.fill( std::size(Set) - 1 );
      Combination end_combination(Set, indexes_end);
      ++end_combination;
      return end_combination;
    }


    

  };
	
	
}//namespace
