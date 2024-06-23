namespace ffd::user_space{

  template<class T1,
	   class T2>
  std::ostream &
  operator << (std::ostream &out,
	       std::pair<T1, T2> const& t){
    return out << std::make_tuple(t.first, t.second);
  }


}//namespace
