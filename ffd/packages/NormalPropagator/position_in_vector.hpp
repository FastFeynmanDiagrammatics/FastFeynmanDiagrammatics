namespace ffd::normal_propagator{

  template<typename T>
  uint
  position_in_vector(T el,
		     std::vector<T> vec){
    for(uint j=0; j<size(vec); ++j){
      if(el == vec[j]){
	return j;
      }
    }
    return size(vec);
  }

}//namespace
