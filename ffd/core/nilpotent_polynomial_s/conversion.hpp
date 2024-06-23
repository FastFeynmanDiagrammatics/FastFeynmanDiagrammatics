namespace ffd::nilpotent_polynomial{

  template <typename Field, SetT SetType>
  NilpotentPolynomial<Field, SetType>::NilpotentPolynomial(ffd::ranked_zeta::RankedZeta<Field, SetType> const & Z){

    auto Z_ = Z;
    int cardinality = Z_.cardinality;
	coef.resize(cardinality);

    for (int j=1; j<cardinality; j*=2){
      for (int S=cardinality-1; S>=j; --S){
        if (j & S){
          BinaryInt alt_S = S-j;
          for (BinaryInt k=0; k<=Z.rank_size; ++k){
            Z_[S+k*cardinality] -=  Z_[alt_S+k*cardinality];
          }
        }
      }
    }

	if (SetType == Full){
		for (int S=0; S<cardinality; ++S){
		  coef[S] = Z_[S+__builtin_popcount(S)*cardinality];
		}
	}
	else if (SetType == Even){
		for (int S=0; S<cardinality; ++S){
		  int pop_S= __builtin_popcount(S);
		  if (pop_S%2==0)
		  {
	  		coef[S] = Z_[S+pop_S/2*cardinality];
		  }
		  else{
			  coef[S] =0.;
		  }
		}
	}

  }

}//namespace
