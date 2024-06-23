namespace ffd::nilpotent_polynomial{

  template<typename Field, SetT SetType>
  std::ostream& operator<<(std::ostream& ret, const NilpotentPolynomial<Field, SetType>& P_){
    if(!P_.NotZero){
      ret<<"ZeroPoly";
    }else if(P_.size() == 1){
      ret<<"scalar("<<P_[0]<<")";
    }else{
      for(BinaryInt V = 0; V < BinaryInt(P_.size()); ++V){
	if(V != 0){
	  ret<<"+(";
	}
	ret<<P_[V];
	if(V != 0){
	  ret<<")";
	}
	for(BinaryInt Sz = V, j=0; Sz!=0; Sz >>= 1, ++j){
	  if(Sz&1){
	    ret<<P_.variable_names[j];
	  }
	}
      }
    }
    return ret;
  }


}//namespace ffd::nilpotent_polynomial
