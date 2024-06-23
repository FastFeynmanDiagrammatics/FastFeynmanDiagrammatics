namespace ffd::user_space{

  template<class coord_d,
	   class field_d>
  auto
  Nambu(QFSum<coord_d, field_d> const& x){
    auto ret = x;
    for(uint j=0; j<size(ret.fields); ++j){
      ret.fields[j].first = Nambu(ret.fields[j].first);
    }
    return ret;
  }


}//namespace
