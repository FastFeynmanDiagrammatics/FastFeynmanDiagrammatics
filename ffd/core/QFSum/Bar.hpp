namespace ffd::user_space{

  template<class coord_d,
	   class field_d>
  auto
  Bar(QFSum<coord_d, field_d> const& x){
    auto ret = x;
    for(uint k=0; k<size(x.plus_pos)-1; ++k){
      auto prod_size = x.plus_pos[k+1]-x.plus_pos[k];
      for(uint j=x.plus_pos[k]; j<x.plus_pos[k+1]; ++j){
	auto j_comp = x.plus_pos[k] + (prod_size-(j-x.plus_pos[k])-1); 
	ret.fields[j] = x.fields[j_comp];
      }
    }
    for(uint j=0; j<size(ret.fields); ++j){
      auto & y = ret.fields[j].first;
      y = Bar(y);
    }
    return ret;
  }

}//namespace
