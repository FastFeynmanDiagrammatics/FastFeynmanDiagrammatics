namespace ffd::normal_principal_blocks{

  template<typename T>
  bool
  is_double_element(int k,
		    std::vector<T> const& v){
    for(uint j=0; j<size(v); ++j){
      if(j != k && v[j] == v[k]) return true; 
    }
    return false;
  }

}//namespace
