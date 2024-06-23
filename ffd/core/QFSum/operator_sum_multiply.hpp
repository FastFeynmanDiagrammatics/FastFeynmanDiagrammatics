namespace ffd::user_space{

  template<class c_1,
  	   class c_2,
  	   class fld1,
  	   class fld2>
  typename std::enable_if_t<std::is_same_v<c_1, nothing_d> ||
  			    std::is_same_v<c_1, c_2>,
  			    QFSum<c_2, decltype(std::declval<fld1>()
					      *std::declval<fld2>())>>
  operator*(QFSum<c_1, fld1> const & q1,
  	    QFSum<c_2, fld2> const & q2){
    QFSum<c_2, decltype(std::declval<fld1>()
			*std::declval<fld2>())> ret;
    ret.plus_pos.resize(size(q1.coef)*size(q2.coef)+1);
    ret.coef.resize(size(q1.coef)*size(q2.coef));
    uint counter_k = 0;
    for(uint k1=0; k1<size(q1.plus_pos)-1; ++k1){
      for(uint k2=0; k2<size(q2.plus_pos)-1; ++k2){
	ret.coef[counter_k] = q1.coef[k1]*q2.coef[k2];
	counter_k++;
	uint counter_j = 0;
	for(uint j1=q1.plus_pos[k1]; j1<q1.plus_pos[k1+1]; ++j1){
	  if constexpr(std::is_same_v<c_1, c_2>){
	      ret.fields.push_back(q1.fields[j1]);
	    }else{
	    ret.fields.push_back(std::make_pair(q1.fields[j1].first,
						q2.fields[0].second));
	  }
	  counter_j++;
	}
	for(uint j2=q2.plus_pos[k2]; j2<q2.plus_pos[k2+1]; ++j2){
	  ret.fields.push_back(q2.fields[j2]);
	  counter_j++;
	}
	ret.plus_pos[counter_k] = counter_j + 
	  ret.plus_pos[counter_k-1];
      }
    }
    return ret;
  }
  
  
  template<typename c_1,
  	   typename fld1,
  	   typename fld2>
  QFSum<c_1, decltype(std::declval<fld1>()
		      *std::declval<fld2>())>
  operator*(fld2 x,
  	    QFSum<c_1, fld1> const & q1){
    QFSum<c_1, decltype(std::declval<fld1>()
			*std::declval<fld2>())> ret;
    ret.fields = q1.fields;
    ret.plus_pos = q1.plus_pos;
    ret.coef.resize(size(q1.coef));
    for(uint j=0; j<size(q1.coef); ++j){
      ret.coef[j] = q1.coef[j]*x;
    }
    return ret;
  }


  template<typename c_1,
  	   typename fld1,
  	   typename fld2>
  QFSum<c_1, decltype(std::declval<fld1>()
		      *std::declval<fld2>())> 
  operator*(QFSum<c_1, fld1> const & q1,
  	    fld2 x){
    return x*q1;
  }


  template<typename c_1,
  	   typename fld1,
  	   typename fld2>
  QFSum<c_1, decltype(std::declval<fld1>()
		      +std::declval<fld2>())>
  operator+(QFSum<c_1, fld1> const & q1,
  	    QFSum<c_1, fld2> const & q2){
    QFSum<c_1, decltype(std::declval<fld1>()
			+std::declval<fld2>())> ret;
    ret.fields = q1.fields;
    for(uint j=0; j<size(q2.fields); ++j){
      ret.fields.push_back(q2.fields[j]);
    }
    ret.plus_pos = q1.plus_pos;
    auto last = ret.plus_pos[size(ret.plus_pos)-1];
    for(uint k=1; k<size(q2.plus_pos); ++k){
      ret.plus_pos.push_back(last+q2.plus_pos[k]);
    }
    ret.coef = q1.coef;
    for(uint j=0; j<size(q2.coef); ++j){
      ret.coef.push_back(q2.coef[j]);
    }
    return ret;
  }

  
  template<typename c_1,
  	   typename fld1>
  QFSum<c_1, fld1>
  operator-(QFSum<c_1, fld1> const & q1){
    return q1*(-1.);
  }


  template<typename c_1,
  	   typename fld1,
  	   typename fld2>
  QFSum<c_1, decltype(std::declval<fld1>()
		      +std::declval<fld2>())>
  operator-(QFSum<c_1, fld1> const & q1,
  	    QFSum<c_1, fld2> const & q2){
    return q1 + (-q2);
  }



  template<typename c_1,
  	   typename fld1,
  	   typename fld2>
  auto
  operator+=(QFSum<c_1, fld1> & q1,
  	     QFSum<c_1, fld2> const & q2){
    return q1 = q1 + q2;
  }
  
  
  template<typename c_1,
  	   typename fld1,
  	   typename fld2>
  auto
  operator-=(QFSum<c_1, fld1> & q1,
  	     QFSum<c_1, fld2> const & q2){
    return q1 = q1 - q2;
  }


  template<typename c_1,
  	   typename fld1,
  	   typename fld2>
  auto
  operator*=(QFSum<c_1, fld1> & q1,
  	     QFSum<c_1, fld2> const & q2){
    return q1 = q1*q2;
  }


}//namespace
