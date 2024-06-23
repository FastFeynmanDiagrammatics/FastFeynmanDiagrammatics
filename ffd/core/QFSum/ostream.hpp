namespace ffd::user_space{

  template<class coord_d,
	   class field_d>
  std::ostream &
  operator << (std::ostream &out,
	       QFSum<coord_d, field_d> const& qfs){
    for(uint k=0; k<size(qfs.plus_pos)-1; ++k){
      out << std::showpos;
      out << qfs.coef[k] << " ";
      out << std::noshowpos;
      for(uint j=qfs.plus_pos[k]; j<qfs.plus_pos[k+1]; ++j){
	auto [qf, X] = qfs.fields[j];
	if(Statistics(qf) == phys::fermi){
	  //	  const char* s = u8"\u0444";
	  out << "\u03A8";
	  //	  out << "psi";
	  //
	  if(Direction(qf) == phys::ou){
	    out << "\u207A";
	    //	out << "\u2020";
	  }
	}
	out << "_{";
	if(Component(qf) == 1){
	  //      out << "\u207
	  out << "\u2191";
	}else if(Component(qf) == -1){
	  out << "\u2193";
	}else{
	  out << std::noshowpos;
	  out << Component(qf);
	  out << std::showpos;
	}
	out << "}";
	if constexpr(!std::is_same_v<coord_d, nothing_d>){
	    out << X << ' ';
	  }
      }
    }
    out << std::noshowpos;
    return out;
  }

}//namespace
