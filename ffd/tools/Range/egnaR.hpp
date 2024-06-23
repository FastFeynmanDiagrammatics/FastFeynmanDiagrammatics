namespace ffd::vector_range{

  template<typename IntType1,
	   typename IntType2 = IntType1>
  std::vector<IntType2>
  egnaR(IntType1 upp, IntType2 low){
    if((long) low >= (long) upp){
      return std::vector<IntType2>();
    }
    std::size_t const n_iter = long(upp) - long(low);
    std::vector<IntType2> ret(n_iter);
    for(std::size_t j = 0; j < n_iter; ++j){
      ret[j] = low +IntType2(long(n_iter)-1-j);
    }
    return ret;
  }


  
  template<typename IntType>
  std::vector<IntType>
  egnaR(IntType upp){
    return egnaR(upp, IntType(0));
  }

  
}//namespace
