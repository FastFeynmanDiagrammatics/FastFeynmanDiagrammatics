namespace ffd::user_space{
  template<int d, typename... args>
  std::ostream & operator << (std::ostream &out,
			      ffd::class_tuple::
			      ClassTupleTypeIntInt<ffd::periodic_coordinate::
			      PeriodicCoordinate, d, 0,
			      args...> const& X)
  {
    out << std::noshowpos;
    out << '(' << ffd::get<0>(X)();
    ///UGLYYYY
    if constexpr(d >= 1 ){
	out << ";" << ffd::get<1>(X)();
      }
    if constexpr(d >= 2 ){
	out << "," << ffd::get<2>(X)();
      }
    if constexpr(d >= 3 ){
	out << "," << ffd::get<3>(X)();
      }
    if constexpr(d >= 4 ){
	out << "," << ffd::get<4>(X)();
      }
    if constexpr(d >= 5 ){
	out << "," << ffd::get<5>(X)();
      }
    if constexpr(sizeof...(args)>d+1 ){
	out << ";" << ffd::get<d+1>(X)();
      }
    
    out << ')';
    return out;
  }


}//namespace
