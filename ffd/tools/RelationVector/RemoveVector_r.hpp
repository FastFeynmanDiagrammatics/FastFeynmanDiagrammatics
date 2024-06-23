namespace ffd::relation_vector{

  template<class T>
  void
  RemoveVector_r(std::size_t block,
		 Container<T>& C){
    std::size_t const size_remove =
      C.split[block+1] - C.split[block];
    for(uint j=C.split[block]; j<size(C.container)-size_remove; ++j){
      C.container[j] = C.container[j+size_remove];
    }
    C.container.resize(size(C.container)-size_remove);
    for(uint j=block+1; j<size(C.split)-1; ++j){
      C.split[j] = C.split[j+1] - size_remove;
    }
    C.split.resize(size(C.split)-1);
  }

}//namespace
