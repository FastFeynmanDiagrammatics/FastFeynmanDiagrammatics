namespace ffd::ranked_zeta{

  template <typename Field, SetT SetType>
  RankedZeta<Field, SetType>::RankedZeta(ffd::nilpotent_polynomial::NilpotentPolynomial<Field, SetType> const & P){

    cardinality = P.size();
  	rank_size = __builtin_popcount(cardinality-1);
    coef.resize(cardinality*(rank_size+1));

    std::vector<BinaryInt> pop(cardinality);
  	for (BinaryInt S=0; S<cardinality; ++S){
  	  pop[S] = __builtin_popcount(S);
  	  coef[S+pop[S]*cardinality] = P[S];
  	}

  	for (int j=1; j<(int)cardinality; j*=2){
  	  for (int S=cardinality-1; S>=j; --S){
  			if (j & S){
  			  BinaryInt alt_S = S-j;
  			  for (BinaryInt k=0; k<=pop[S]; ++k){
  				  coef[S+k*cardinality] += coef[alt_S+k*cardinality];
  			  }
  			}
  	  }
  	}

    if (SetType==Even){
      rank_size/=2;
      std::vector<Field> new_coef (cardinality*(rank_size+1),0.);
      for (int S=cardinality-1; S>=0; --S){
        for (int k=0; k<=(int)pop[S]; k+=2){
          new_coef[S+k*cardinality/2] = coef[S+k*cardinality];
        }
      }
      coef = new_coef;
    }

  }

  // template <typename Field>
  // RankedZeta<Field, Even>::RankedZeta(RankedZeta<Field, Full> const & Z){
	// cardinality = Z.cardinality;
	// rank_size = Z.ranksize/2;
	// coef.resize(cardinality*(rank_size+1));
	// for (int S=cardinality-1; S>=0; --S){
	//   for (int k=0; k<= __builtin_popcount(S); k+=1){
	// 	if(k%2==1){
	// 		assert(Z[S+k*cardinality]==0);
	// 	}
	// 	else{
	// 		coef[S+k*cardinality/2] = Z[S+k*cardinality];
	// 	}
	//   }
	// }
  // }
  //
  // template <typename Field>
  // RankedZeta<Field, Full>::RankedZeta(RankedZeta<Field, Even> const & Z){
  // cardinality = Z.cardinality;
  // rank_size = Z.ranksize*2;
  // coef.resize(cardinality*(rank_size+1));
  // for (int S=cardinality-1; S>=0; --S){
	// for (int k=0; k<= __builtin_popcount(S); k+=2){
	//   coef[S+k*cardinality] = Z[S+k*cardinality/2];
	// }
  // }
  // }

}//namespace
