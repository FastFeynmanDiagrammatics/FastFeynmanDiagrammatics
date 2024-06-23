namespace ffd::vector_range{

  template<typename IntType2,
	   typename IntType1 = IntType2>
  std::vector<IntType1>
  Range(IntType1 low, IntType2 upp){
    if((long) low >= (long) upp){
      return std::vector<IntType1>();
    }
    long low_l = low;
    std::size_t const size = long(upp) - low_l;
    std::vector<IntType1> ret(size);
    for(std::size_t j=0; j<size; ++j){
      ret[j] = low+j;
    }
    return ret;
  }

  
  template<typename IntType>
  std::vector<IntType>
  Range(IntType upp){
    return Range(IntType(0), upp);
  }
  
	
}//namespace
