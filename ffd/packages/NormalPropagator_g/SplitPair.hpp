namespace ffd::user_space{

  template<class coord_d, class field_d>
  auto
  SplitPair(QFSum<coord_d, field_d> const& x){
    assert((size(x.fields) == 2));
    return std::make_pair(x.fields[0], x.fields[1]);
  }

}//namespace
