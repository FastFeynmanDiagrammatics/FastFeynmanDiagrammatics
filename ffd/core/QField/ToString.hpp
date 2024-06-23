namespace ffd::user_space{

  auto
  ToString(QField qf){
    std::stringstream out;
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
    return out.str();
  }

}//namespace
