namespace ffd{
  template<int j, typename T>
  T& get(ffd::class_tuple::ShellStruct<j, T>& t_){
    return t_.Content;
  }                            

  template<int j, typename T>
  T const& get(ffd::class_tuple::ShellStruct<j, T> const& t_){
    return t_.Content;
  }                            

  // template<int j, template<typename> class C, typename T>
  // const C<T>& Get(const ffd::class_tuple::ShellStruct<j, C<T>>& t_){
  //   return t_.Content;
  // }
  
}//namespace ffd

