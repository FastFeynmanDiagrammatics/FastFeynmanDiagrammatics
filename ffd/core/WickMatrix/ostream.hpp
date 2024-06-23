namespace ffd::user_space{

  template<typename T>
  std::ostream&
  operator<<(std::ostream& out, const ffd::wick_matrix::WickMatrix<T>& M_){
    out<<SetOstreamFloatFlags();
    for(int j=0; j < size(M_); ++j){
      out<<"|";
      for(int k=0; k < size(M_); ++k){
	out<<M_(j, k);
	if(k != size(M_)-1){
	  out<<" ";
	}
      }
      out<<"|\n";
    }
    return out;
  }

}//namespace ffd::user_space
