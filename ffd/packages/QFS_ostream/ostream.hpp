namespace ffd::user_space{

  template<typename coord_t, typename fld>
  std::ostream & operator << (std::ostream &out,
			      QFS<coord_t, fld> const& q)
  {
    out << std::showpos;
    for(uint j=0; j<size(q); ++j){
      // if(std::is_same_v<Complex, fld> && j!=0){
      // 	out<<"+";
      // }
      out << std::showpos;
      out << q.coef[j] << " ";
      for(uint k=0; k<size(q[j]); ++k){
	auto [qf, X] = q[j][k];
	out << qf ;
	if constexpr(!std::is_same_v<coord_t, void_t>){
	    out << X << ' ';
	  }
      }
      std::cerr << '\n';
    }
    out << std::noshowpos;
    return out;
  }
  
}//namespace
