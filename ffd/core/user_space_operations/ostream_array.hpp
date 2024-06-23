namespace ffd::user_space{

  template<class T, std::size_t d>
  std::ostream &
  operator << (std::ostream &out,
	       std::array<T, d> const& a){
    out << "(";
    for(uint j=0; j<size(a)-1; ++j){
      out << a[j] << ", ";
    }
    if(size(a) > 1){
      out << a[size(a)-1];
    }
    out << ")";
    return out;
  }
  
}//namespace
