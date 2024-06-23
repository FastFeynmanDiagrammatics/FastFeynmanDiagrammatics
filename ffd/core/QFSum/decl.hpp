namespace ffd::user_space{
  struct nothing_d{};

  template<class coord_d = nothing_d,
	   class field_d = Real>
  struct QFSum{
    std::vector<std::pair<QField, coord_d>> fields;
    std::vector<std::size_t> plus_pos;
    std::vector<field_d> coef;

    QFSum(){plus_pos.push_back(0ul);plus_pos.push_back(0ul);coef.push_back(1.);}

    QFSum
    operator()(int comp) const{
      auto ret = *this;
      for(std::size_t j = plus_pos[0];
	  j<plus_pos[size(plus_pos)-1];
	  ++j){
	fields[j].first =
	  SetComponent(fields[j].first, comp);
      }
      return ret;
    }
        
    QFSum
    operator()(char name_component) const{
      if( name_component == 'u')
	return (*this)(1);
      if( name_component == 'd')
	return (*this)(-1);
      return (*this)(0);
    }
    
    template<typename coord2_d>
    typename std::enable_if_t<std::is_same_v<coord_d, nothing_d>,
			      QFSum<coord2_d, field_d>>
    operator()(coord2_d const&) const;
    
  };
  
}//namespace
