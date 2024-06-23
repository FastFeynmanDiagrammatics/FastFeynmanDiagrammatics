namespace ffd::itime_lattice_proposer{

  template<class... T>
  auto
  split_arguments(T&&... arg){
    return std::make_tuple(std::forward<T>(arg)...);
  }

}//namespace
