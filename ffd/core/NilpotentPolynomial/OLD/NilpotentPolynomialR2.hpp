#ifndef FFD_POLY2_HEADER_NOT_BEEN_HERE
#define FFD_POLY2_HEADER_NOT_BEEN_HERE


namespace ffd::nilpotent_polynomial{
  
  template <int degree, typename Field = Real>
  struct NilpotentPolynomialR2{

    enum class PolynomialState{Zero, Scalar, Polynomial};

    PolynomialState State = PolynomialState::Polynomial;
    
    //the number of coefficients of the polynomials
    static constexpr BinaryInt NumberCoef = ffd::core_math::Pow1ToInt2(2, degree);

    //this is used just for printing
    static const char variable_names[26];
    
    //the subset of \{\xi_1, ..., \xi_n\} in which we perform operations
    BinaryInt set = NumberCoef-1;
    
    //the array of coefficients
    Field coef[NumberCoef] = {0.};

    //static BinaryInt SetShiftsArray[degree];
    
    NilpotentPolynomialR2() {}
    
    NilpotentPolynomialR2(const char* WhichPolynomial){
      //std::cerr<<"STRN";
      for(BinaryInt V=0; V<NumberCoef; ++V){
	coef[V] = 0;
      }
      if(WhichPolynomial == "Zero"){
	State = PolynomialState::Zero;
      }else if(WhichPolynomial == "One"){
	State = PolynomialState::Scalar;
	coef[0] = 1.;
      }else{
	throw std::invalid_argument("NilpotentPolynomialR2::NilpotentPolynomialR2(std::string): argument unknown");
      }
    }

    
    NilpotentPolynomialR2(Field x){
      for(BinaryInt V=1; V<NumberCoef; ++V){
	coef[V] = 0.;
      }
      coef[0] = x;
      State = PolynomialState::Scalar;
    }
    
    
    NilpotentPolynomialR2& operator=(const NilpotentPolynomialR2& source){
      for(BinaryInt V=0; V<NumberCoef; ++V){
    	coef[V] = source.coef[V];
      }
      set = source.set;
      State = source.State;
      //std::cerr<<"EQUL";
      return *this;
    }
    
    
    NilpotentPolynomialR2(const NilpotentPolynomialR2& source){
      //std::cerr<<"COPY";
      *this = source;
    }
    
  
    //returns the maximal value of coef
    Real MaxAbsCoef(){
      Real abs_err_max = 0;
      for(BinaryInt V=0; V<NumberCoef; ++V){
	if(abs_err_max < std::abs(coef[V])){
	  abs_err_max = std::abs(coef[V]);
	}
      }
      return abs_err_max;
    }
  
  
  
    //it converts a set S belonging to the internal convention
    //(which means that S \in {0, 1, 2^degree-1})
    //in a set belonging to the external convention,
    //which is determined by "set"
    inline BinaryInt ConvertInternalSetToExternal(BinaryInt S) const{
      BinaryInt set_copy = set, V = 0;
      for(int j=0; set_copy>0; ++j){
	if(set_copy%2 > 0){
	  V += (S%2) << j;
	  S >>= 1;
	}
	set_copy >>= 1;
      }
      return V;
    }

  
  
    //it converts a set belonging to the external convention
    //to a set belonging to the internal convention
    //for coefficients
    BinaryInt ConvertExternalSetToInternal(BinaryInt S){
      BinaryInt set_copy = set, V = 0;
      for(BinaryInt pow_NumericBase = 1; S>0; S /= 2){
	V += (S%2)*pow_NumericBase;
	if(set_copy%2 > 0){
	  pow_NumericBase *= 2;
	}
	set_copy/= 2;
      }
      return V;
    }
  

    //the () operator uses the external representation
    Field& operator()(BinaryInt S){
      return coef[ConvertExternalSetToInternal(S)];
    }

    
    //the [] operator uses the internal representation
    inline Field& operator[](const BinaryInt& S){ return coef[S];};


    //the [] operator uses the internal representation
    inline const Field& operator[](const BinaryInt& S) const{ return coef[S];};

    NilpotentPolynomialR2& operator+=(Field number){
      if(State == PolynomialState::Zero){
	State = PolynomialState::Scalar;
      }
      coef[0] += number;
      return *this;
    }

    NilpotentPolynomialR2& operator-=(Field number){
      return *this += -number;
    }
    
    NilpotentPolynomialR2 operator+(Field number) const{
      NilpotentPolynomialR2 ret = *this;
      return ret += number;
    }

    
    NilpotentPolynomialR2 operator-(Field number) const{
      NilpotentPolynomialR2 ret = *this;
      return ret -= number;
    }

    
    NilpotentPolynomialR2& operator*=(Field number){
      if(State != PolynomialState::Zero){
	for(BinaryInt S=0; S<NumberCoef; S++){
	  coef[S] *= number;
	}
      }
      return *this;
    }

    
    NilpotentPolynomialR2 operator*(Field number) const{
      NilpotentPolynomialR2 ret = *this;
      return ret *= number;
    }

    
    //with this we can write -P[[\xi]]
    NilpotentPolynomialR2 operator-() const{
      NilpotentPolynomialR2 ret = *this;
      return ret *= -1.;
    }

    NilpotentPolynomialR2& operator/=(Field number){
      *this *= 1./number;
      return *this;
    }

    NilpotentPolynomialR2 operator/(Field number) const{
      NilpotentPolynomialR2 ret = *this;
      return ret /= number;
    }
    
    
    //The vector product multiply each \xi_j by some vector component v_j  
    NilpotentPolynomialR2 operator*(std::vector<Field> weight) const{
      NilpotentPolynomialR2 ret = *this;
      for(BinaryInt S=0; S<NumberCoef; S++){
	for(int j=0; j<degree; ++j){
	  BinaryInt pow_NumericBase = ffd::core_math::Pow1ToInt2(2, j);
	  for(int l=1; l<2; ++l){
	    if(ffd::set_theory::IsSet1SubsetOfSet2(pow_NumericBase*l, ConvertInternalSetToExternal(S))){
	      ret[S] *= weight[j];
	    }
	  }
	}
      }
    }

    
    NilpotentPolynomialR2& operator*=(std::vector<Field> weight){
      return *this = *this * weight;
    }

    
    NilpotentPolynomialR2& operator+=(const NilpotentPolynomialR2& add){
      if(add.State == PolynomialState::Zero){
	return *this;
      }
      if(add.State == PolynomialState::Scalar){
	return *this += coef[0];
      }
      State = PolynomialState::Polynomial;
      set = add.set;
      for(BinaryInt S=0; S<NumberCoef; ++S){
	coef[S] += add[S];
      }
      return *this;
    }


    NilpotentPolynomialR2& operator-=(const NilpotentPolynomialR2& add){
      return *this += -add;
    }

    
    //addition between polynomials
    NilpotentPolynomialR2 operator+(const NilpotentPolynomialR2& A1) const{
      NilpotentPolynomialR2 ret = *this;
      return ret += A1;
    }

    NilpotentPolynomialR2 operator-(const NilpotentPolynomialR2& A1) const{
      NilpotentPolynomialR2 ret = -A1;
      return ret += *this;
    }


    // inline BinaryInt SetShifts(BinaryInt V){
    //   BinaryInt ret = 0;
    //   for(BinaryInt j=0; V > 0; ++j, V >>= 1){
    // 	SetShiftsArray[ret] = j;
    // 	ret += V&1;
    //   }
    //   return ret;
    // }

    
    
    
    //multiply *this by a polynomial 
    NilpotentPolynomialR2 operator*(const NilpotentPolynomialR2& fac_) const{
      using ffd::set_theory::BinaryShiftsOf1In2;
      if(State == PolynomialState::Zero || fac_.State == PolynomialState::Zero){
	return "Zero";
      }else if(State == PolynomialState::Scalar){
	return fac_*coef[0];
      }else if(fac_.State == PolynomialState::Scalar){
	return *this * fac_[0];
      }
      NilpotentPolynomialR2 ret = *this * fac_[0] + fac_ * coef[0];

      BinaryInt Subsets[NumberCoef];

      for(BinaryInt V=1; V<NumberCoef; ++V){
	//	const BinaryInt card_V = ret.SetShifts(V);
	//	const BinaryInt2 card_V = BinaryShiftsOf1In2<degree>(V, ArrayTemp);
	BinaryInt NumSubsets = 0;
	for(BinaryInt S = 1; S < V; ++S){
	  Subsets[NumSubsets] = S;
	  NumSubsets += ((S&(~V))==0);
	}
	for(BinaryInt S = 0; S < NumSubsets; ++S){
	  ret[V] += coef[Subsets[S]] * fac_[V-Subsets[S]];
	}
      }
      return ret;
    }
    
    
    //multiplication between polynomials
    NilpotentPolynomialR2& operator*=(const NilpotentPolynomialR2& fac){
      return *this = *this * fac;
    }
    
    // //divide *this by a polynomial 
    // NilpotentPolynomialR2 operator/(const NilpotentPolynomialR2& den) const{
    //   using ffd::set_theory::BinaryShiftsOf1In2;
    //   if(State == PolynomialState::Zero){
    // 	return "Zero";
    //   }
    //   NilpotentPolynomialR2 ret = *this / den[0];
    //   if(den.State == PolynomialState::Scalar){
    // 	return ret;
    //   }
    //   NilpotentPolynomialR2 den_norm = den / den[0];
    //   BinaryInt2 ArrayTemp[degree];
    //   for(BinaryInt V=1; V<NumberCoef; ++V){
    // 	//	const BinaryInt card_V = ret.SetShifts(V);
    // 	const BinaryInt2 card_V = BinaryShiftsOf1In2<degree>(V, ArrayTemp);
    // 	for(BinaryInt S = 0, S_V = 0, S_copy = 0; S < (1<<card_V)-1; ++S, S_copy=S, S_V=0 ){
    // 	  //	  BinaryInt S_V = 0, S_copy = S;//, S_V_c=0, V_S = V-S;
    // 	  for(BinaryInt2 k=0; k < card_V; ++k, S_copy >>= 1){//, V_S >>= 1){
    // 	    //	    S_V += (S_copy&1)<<ret.SetShiftsArray[k];
    // 	    S_V += (S_copy&1)<<ArrayTemp[k];
    // 	    //	    S_V_c += (V_S%2)<<ret.SetShiftsArray[k];
    // 	  }
    // 	  //	  std::cerr<<S_V<<std::endl;
    // 	  ret[V] -= ret[S_V]*den_norm[V-S_V];
    // 	}
    //   }
    //   return ret;
    // }
    

    NilpotentPolynomialR2& operator/=(const NilpotentPolynomialR2& den){
      return *this = *this / den;
    }
    
    
    //restrict the polynomial over a specified set S (that can be of lower order r)
    //and it stores it in the second argument
    //WARNING: it works only if nilpotency == 1
    template<int r>
    void RestrictTo1StoreIn2(BinaryInt S, NilpotentPolynomialR2<r, Field>& dest){
      dest.State = State;
      dest = 0;
      dest.set = S;
      if(set==0){
	dest.set = 0;
	dest[0] = coef[0];
      }else if(S==0){
	dest[0] = coef[0];
      }else{
	for(BinaryInt V=0; V < (1<<r); ++V){
	  dest[V] = coef[  ConvertExternalSetToInternal( dest.ConvertInternalSetToExternal(V) )  ];
	}
      }
    }
  };



  template <int degree, typename Field>
  const char NilpotentPolynomialR2<degree, Field>::variable_names[] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'};


  // template<int degree, typename Field>
  // BinaryInt NilpotentPolynomialR2<degree, Field>::SetShiftsArray[] = {0};

  
}//namespace ffd::nilpotent_polynomial




#endif
