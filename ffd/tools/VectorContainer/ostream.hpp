namespace ffd::user_space{

  template<class T>
  std::ostream &
  operator << (std::ostream &out,
	       VectorContainer<T> const& C){
    out << "|";
    for(uint k=0; k<size(C.split)-1; ++k){
      out << ' ';
      for(uint j=C.split[k]; j<C.split[k+1]; ++j){
	out << C.container[j] << ' '; 
      }
      out << '|';
    }
    return out;
  }


}//namespace
