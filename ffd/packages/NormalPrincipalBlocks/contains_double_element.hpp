namespace ffd::normal_principal_blocks{

  template<typename T>
  bool
  contains_double_element(std::vector<T> const& v){
    for(uint j=0; j<size(v); ++j){
      for(uint k=j+1; k<size(v); ++k){
	if(v[j] == v[k]) return true;
      }
    }
    return false;
  }

}//namespace
