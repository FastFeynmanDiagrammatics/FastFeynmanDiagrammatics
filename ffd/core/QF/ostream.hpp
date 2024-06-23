namespace ffd::user_space{

  std::ostream & operator << (std::ostream &out,
			      QF const& qf)
  {
    if(qf.statistics == phys::fermi){
      //	  const char* s = u8"\u0444";
      out << "\u03A8";;
	  //	  out << "psi";

	     //
      if(qf.direction == phys::ou){
	out << "\u207A";
	//	out << "\u2020";
      }
    }
    out << "_{";
    if(qf.component == 1){
      //      out << "\u207
      out << "\u2191";
    }else if(qf.component == -1){
      out << "\u2193";
    }else{
      out << std::noshowpos;
      out << qf.component;
      out << std::showpos;
	}
    out << "}";
    return out;
  }
  


}//namespace
