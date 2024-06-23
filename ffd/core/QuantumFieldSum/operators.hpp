namespace ffd::user_space{

  template<typename coord_t, typename field>
  template<typename coord_t2> 
  typename std::enable_if_t<
			    std::is_same_v<coord_t, void_t>,
			    QuantumFieldSum<coord_t2, field>>
  QuantumFieldSum<coord_t, field>::operator()(coord_t2 const& Y) const{
    QuantumFieldSum<coord_t2, field> ret;
    for(uint j=0; j<size(*this); ++j){
      std::vector<std::pair<QF, coord_t2>> qf_prod;
      for(uint k=0; k<size((*this)[j]); ++k){
	auto qf_now = (*this)[j][k].first;
	qf_prod.push_back(std::make_pair(qf_now, Y));
      }
      ret.fields.push_back(qf_prod);
    }
    ret.coef = this->coef;
    return ret;
  }

    
  template<typename c_1,
  	   typename c_2,
  	   typename fld1,
  	   typename fld2>
  typename std::enable_if_t<std::is_same_v<c_1, void_t> ||
  			    std::is_same_v<c_1, c_2>,
  			    QFS<c_2, decltype(std::declval<fld1>()
					      *std::declval<fld2>())>>
  operator*(QFS<c_1, fld1> const & q1,
  	    QFS<c_2, fld2> const & q2){
    QFS<c_2, decltype(std::declval<fld1>()
		      *std::declval<fld2>())> ret;
    for(uint j1=0; j1<size(q1); ++j1){
      for(uint j2=0; j2<size(q2); ++j2){
	ret.coef.push_back(q1.coef[j1]*q2.coef[j2]);
	std::vector<std::pair<QF, c_2>> q_merged;
	for(uint j3=0; j3<size(q1[j1]); ++j3){
	  auto x = q1[j1][j3];
	  std::pair<QF, c_2> xx;
	  xx.first = x.first;
	  if constexpr(std::is_same_v<void_t, c_1>){
 xx.second = q2[j2][0].second;
}else{
	    xx.second = x.second;
	  }
	  q_merged.push_back(xx);
	}
	for(uint j3=0; j3<size(q2[j2]); ++j3){
	  q_merged.push_back(q2[j2][j3]);
	}
	ret.fields.push_back(q_merged);
      }
    }
    return ret;
  }
  
  
  template<typename c_1,
  	   typename fld1,
  	   typename fld2>
  QFS<c_1, decltype(std::declval<fld1>()
  		    *std::declval<fld2>())>
  operator*(fld2 x,
  	    QFS<c_1, fld1> const & q1){
    QFS<c_1, decltype(std::declval<fld1>()
		      *std::declval<fld2>())> ret;
    ret.fields = q1.fields;
    for(uint j=0; j<size(q1); ++j){
      ret.coef.push_back(q1.coef[j]*x);
    }
    return ret;
  }


  template<typename c_1,
  	   typename fld1,
  	   typename fld2>
  QFS<c_1, decltype(std::declval<fld1>()
  		    *std::declval<fld2>())> 
  operator*(QFS<c_1, fld1> const & q1,
  	    fld2 x){
    return x*q1;
  }


  template<typename c_1,
  	   typename fld1,
	   typename fld2>
  QFS<c_1, decltype(std::declval<fld1>()
		    +std::declval<fld2>())>
  operator+(QFS<c_1, fld1> const & q1,
  	    QFS<c_1, fld2> const & q2){
    QFS<c_1, decltype(std::declval<fld1>()
		      +std::declval<fld2>())> ret;
    ret.fields = q1.fields;
    for(uint j=0; j<size(q2); ++j){
      ret.fields.push_back(q2[j]);
    }

    for(uint j=0; j<size(q1); ++j){
      ret.coef.push_back(q1.coef[j]);
    }
    for(uint j=0; j<size(q2); ++j){
      ret.coef.push_back(q2.coef[j]);
    }
    return ret;
  }

  
  template<typename c_1,
  	   typename fld1>
  QFS<c_1, fld1>
  operator-(QFS<c_1, fld1> const & q1){
    return q1*(-1.);
  }


  template<typename c_1,
  	   typename fld1,
	   typename fld2>
  QFS<c_1, decltype(std::declval<fld1>()
		    +std::declval<fld2>())>
  operator-(QFS<c_1, fld1> const & q1,
  	    QFS<c_1, fld2> const & q2){
    return q1 + (-q2);
  }



  template<typename c_1,
  	   typename fld1,
	   typename fld2>
  auto
  operator+=(QFS<c_1, fld1> & q1,
	     QFS<c_1, fld2> const & q2){
    return q1 = q1 + q2;
  }


  template<typename c_1,
  	   typename fld1,
	   typename fld2>
  auto
  operator-=(QFS<c_1, fld1> & q1,
	     QFS<c_1, fld2> const & q2){
    return q1 = q1 - q2;
  }


  template<typename c_1,
  	   typename fld1,
	   typename fld2>
  auto
  operator*=(QFS<c_1, fld1> & q1,
	     QFS<c_1, fld2> const & q2){
    return q1 = q1*q2;
  }

}//namespace
