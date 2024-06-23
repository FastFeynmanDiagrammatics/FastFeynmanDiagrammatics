namespace ffd::normal_principal_blocks{

  template<typename T>
  void
  remove_double_elements(std::vector<T>& v){
    std::sort(v.begin(), v.end()); 
    auto last = std::unique(v.begin(), v.end());
    v.erase(last, v.end());
  }
  
}//namespace
