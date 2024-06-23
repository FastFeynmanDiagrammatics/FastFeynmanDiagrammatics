namespace ffd::enumerate{

  template<typename vector_t>

  auto
  Enumerate(vector_t vector_v){
    std::vector<std::pair<std::size_t,
			  typename std::decay<decltype(vector_v[0])>::type>> ret( size(vector_v) );

    
    for( auto j: ffd::vector_range::RangeSize( vector_v ) ){
      ret[j] = std::make_pair(std::size_t(j), vector_v[j]);
    }


    return ret;
  }
  
}//namespace
