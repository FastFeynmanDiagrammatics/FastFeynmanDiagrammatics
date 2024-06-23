namespace ffd::user_space{

  template<typename coord_t,
	   typename field>
  QFS<coord_t, field>
  Bar(QFS<coord_t, field> const& qfs){
    QFS<coord_t, field> ret = qfs;

    
    for(uint j=0; j<size(qfs); ++j){
      std::reverse(std::begin(ret[j]),
		   std::end(ret[j]));
    }
    for(uint j=0; j<size(qfs); ++j){
      for(uint k=0; k<size(qfs[j]); ++k){
	ret[j][k].first = Bar(ret[j][k].first);
      }
    }
    
    return ret;
  }
  
}//namespace
