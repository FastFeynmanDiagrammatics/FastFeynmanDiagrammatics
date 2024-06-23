namespace ffd::cartesian_product{

  
  template<auto d,
	   template<typename, typename...> typename SetType,
	   typename ValueType,
	   typename... Args>

  struct CartesianProduct{

    std::array<SetType<ValueType, Args...>, d> Sets;
    
    std::array<std::size_t, d> indexes;
    
    
    
    explicit CartesianProduct(std::array<SetType<ValueType, Args...>, d> Sets_,
			      std::array<std::size_t, d> indexes_):
      Sets(Sets_), indexes(indexes_) {}
    

    explicit CartesianProduct(std::array<SetType<ValueType, Args...>, d> Sets_):
      Sets(Sets_) { indexes.fill(0ul); }


    
    CartesianProduct& operator++(){
      for( auto j: ffd::vector_range::egnaR(d) ){
    	if(  j != 0  &&
    	     indexes[j] == size(Sets[j]) - 1  ){
    	  indexes[j] = 0;
    	}else{
    	  indexes[j] = indexes[j] + 1;
	  break;
    	}
      }
      return *this;
    }


    
    bool operator!=(CartesianProduct const& it) const{
      return indexes != it.indexes;
    }

    
    std::array<ValueType, d>
    operator*() const{
      std::array<ValueType, d> Values;
      for( int j: ffd::vector_range::Range(d) ){
	auto iterator_set = Sets[j].cbegin();
	std::advance(iterator_set, indexes[j]);
	Values[j] = *iterator_set;
      }
      return Values;
    }
    

    CartesianProduct begin() const{
      return CartesianProduct(Sets);
    }

    
    CartesianProduct end() const{
      using ffd::vector_range::Range;
      std::array<std::size_t, d> indexes_end;
      for( auto j: Range(d) ){
	indexes_end[j] = size(Sets[j]) - 1;
      }
      
      CartesianProduct cartesian_end(Sets, indexes_end);
      ++cartesian_end;
      return cartesian_end;
    }

    
  };
	
	
}//namespace
