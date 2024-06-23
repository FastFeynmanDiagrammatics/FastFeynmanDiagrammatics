namespace ffd::user_space{

  template<class coord_d, class field_d>
  template<class coord2_d> 
  typename std::enable_if_t<std::is_same_v<coord_d, nothing_d>,
			    QFSum<coord2_d, field_d>>
  QFSum<coord_d, field_d>::operator()(coord2_d const& Y) const{
    QFSum<coord2_d, field_d> ret;
    ret.coef = this->coef;
    ret.plus_pos = this->plus_pos;
    for(uint j=0; j<size(this->fields); ++j){
      ret.fields.push_back(std::make_pair(this->fields[j].first,
					  Y));
    }
    return ret;
  }
  
}//namespace
