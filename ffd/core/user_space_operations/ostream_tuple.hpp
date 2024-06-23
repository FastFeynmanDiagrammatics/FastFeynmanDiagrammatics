namespace ffd::user_space{

  template<int r,
	   class... T>
  std::ostream &
  ostream_rec_tuple (std::ostream &out,
		     std::tuple<T...> const& t){
    if constexpr(r==0){
	out << "(";
      }
    if constexpr(r < sizeof...(T) ){
	out << std::get<r>(t);
	if constexpr(r+1 < sizeof...(T)){
	  out << "; ";
	}
	return ostream_rec_tuple<r+1>(out, t);
      }else{
      out << ")";
      return out;
    }
  }


  template<class... T>
  std::ostream &
  operator << (std::ostream &out,
	       std::tuple<T...> const& t){
    return ostream_rec_tuple<0>(out, t);
  }
}//namespace
