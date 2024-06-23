

namespace ffd::nilpotent_polynomial{
  
  template <int degree, typename Field = Real>
  class NilpotentPolynomialF{
  public:
    
    //the number of coefficients of the polynomials
    static constexpr BinaryInt NumberCoef = ffd::core_math::Pow1ToInt2(2, degree);

    //this is used just for printing
    static const char variable_names[26];
    
    //the subset of \{\xi_1, ..., \xi_n\} in which we perform operations
    BinaryInt set = NumberCoef-1;
    
    //the array of coefficients
    Field coef[NumberCoef] = {0.};

    enum class PolynomialState{Zero, Scalar, Polynomial};
    PolynomialState State = PolynomialState::Polynomial;
    

    //static BinaryInt SetShiftsArray[degree];
    
    NilpotentPolynomialF() {}
    
    NilpotentPolynomialF(const char* WhichPolynomial){
      //std::cerr<<"STRN";
      for(BinaryInt V=0; V<NumberCoef; ++V){
	coef[V] = 0;
      }
      if(WhichPolynomial == std::string("Zero")){
	State = PolynomialState::Zero;
      }else if(WhichPolynomial == std::string("One")){
	State = PolynomialState::Scalar;
	coef[0] = 1.;
      }else{
	throw std::invalid_argument("NilpotentPolynomialF::NilpotentPolynomialF(std::string): argument unknown");
      }
    }

    
    NilpotentPolynomialF(Field x){
      for(BinaryInt V=1; V<NumberCoef; ++V){
	coef[V] = 0.;
      }
      coef[0] = x;
      State = PolynomialState::Scalar;
    }
    
    
    NilpotentPolynomialF& operator=(const NilpotentPolynomialF& source){
      for(BinaryInt V=0; V<NumberCoef; ++V){
    	coef[V] = source[V];
      }
      set = source.set;
      State = source.State;
      //std::cerr<<"EQUL";
      return *this;
    }
    
    
    NilpotentPolynomialF(const NilpotentPolynomialF& source){
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

    NilpotentPolynomialF& operator+=(Field number){
      if(State == PolynomialState::Zero){
	State = PolynomialState::Scalar;
      }
      coef[0] += number;
      return *this;
    }

    NilpotentPolynomialF& operator-=(Field number){
      return *this += -number;
    }
    
    NilpotentPolynomialF operator+(Field number) const{
      NilpotentPolynomialF ret = *this;
      return ret += number;
    }

    
    NilpotentPolynomialF operator-(Field number) const{
      NilpotentPolynomialF ret = *this;
      return ret -= number;
    }

    
    NilpotentPolynomialF& operator*=(Field number){
      if(State != PolynomialState::Zero){
	for(BinaryInt S=0; S<NumberCoef; S++){
	  coef[S] *= number;
	}
      }
      return *this;
    }

    
    NilpotentPolynomialF operator*(Field number) const{
      NilpotentPolynomialF ret = *this;
      return ret *= number;
    }

    
    //with this we can write -P[[\xi]]
    NilpotentPolynomialF operator-() const{
      NilpotentPolynomialF ret = *this;
      return ret *= -1.;
    }

    NilpotentPolynomialF& operator/=(Field number){
      *this *= 1./number;
      return *this;
    }

    NilpotentPolynomialF operator/(Field number) const{
      NilpotentPolynomialF ret = *this;
      return ret /= number;
    }
    
    
    //The vector product multiply each \xi_j by some vector component v_j  
    NilpotentPolynomialF operator*(std::vector<Field> weight) const{
      NilpotentPolynomialF ret = *this;
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

    
    NilpotentPolynomialF& operator*=(std::vector<Field> weight){
      return *this = *this * weight;
    }

    
    NilpotentPolynomialF& operator+=(const NilpotentPolynomialF& add){
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


    NilpotentPolynomialF& operator-=(const NilpotentPolynomialF& add){
      return *this += -add;
    }

    
    //addition between polynomials
    NilpotentPolynomialF operator+(const NilpotentPolynomialF& A1) const{
      NilpotentPolynomialF ret = *this;
      return ret += A1;
    }

    NilpotentPolynomialF operator-(const NilpotentPolynomialF& A1) const{
      NilpotentPolynomialF ret = -A1;
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
    NilpotentPolynomialF operator*(const NilpotentPolynomialF& fac_) const{
      if(State == PolynomialState::Zero || fac_.State == PolynomialState::Zero){
	return "Zero";
      }else if(State == PolynomialState::Scalar){
	return fac_ * coef[0];
      }else if(fac_.State == PolynomialState::Scalar){
	return *this * fac_[0];
      }
      
      NilpotentPolynomialF ret__ = *this*fac_[0];
      for(BinaryInt V=1; V<NumberCoef; ++V){
	for(BinaryInt S = V; S != 0; S = ((S-1)&V)){
	  ret__[V] += fac_[S] * coef[V-S];
	}
      }
      return ret__;
    }
    
    
    //multiplication between polynomials
    NilpotentPolynomialF& operator*=(const NilpotentPolynomialF& fac_){
      return *this = *this * fac_;
    }


    //divide *this by a polynomial 
    NilpotentPolynomialF operator/(const NilpotentPolynomialF& den_) const{
      if(State == PolynomialState::Zero){
	return "Zero";
      }else if(den_.State == PolynomialState::Scalar){
	return *this/den_[0];
      }else if(den_.State == PolynomialState::Zero){
	throw std::invalid_argument("NilpotentPolynomialF::operator/: division by zero");
      }

      NilpotentPolynomialF ret__ = *this/den_[0];
      NilpotentPolynomialF den__ = den_/den_[0];
      for(BinaryInt V = 1; V < NumberCoef; ++V){
	for(BinaryInt S = V; S != 0; S = ((S-1)&V)){
	  ret__[V] -= den__[S] * ret__[V-S];
	}
      }
      return ret__;
    }


    NilpotentPolynomialF& operator/=(const NilpotentPolynomialF& den){
      return *this = *this / den;
    }
    
    
    //restrict the polynomial over a specified set S (that can be of lower order r)
    //and it stores it in the second argument
    //WARNING: it works only if nilpotency == 1
    template<int r>
    void RestrictTo1StoreIn2(BinaryInt S, NilpotentPolynomialF<r, Field>& dest){
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
  const char NilpotentPolynomialF<degree, Field>::variable_names[] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'};


  // template<int degree, typename Field>
  // BinaryInt NilpotentPolynomialF<degree, Field>::SetShiftsArray[] = {0};

  template<BinaryInt S, typename Field>
  Field ShiftByXi(Field x_){
    return x_;
  }


  template<BinaryInt V, int degree, typename Field>
  NilpotentPolynomialF<degree-ffd::set_theory::CardinalitySet(V), Field>
  ShiftByXi(const NilpotentPolynomialF<degree, Field>& P_){
    if(P_.set != P_.NumberCoef - 1){
      throw std::invalid_argument("ffd::nilpotent_polynomial::ShiftByXi: set of argument must be maximal");
    }
    constexpr BinaryInt set = (1<<degree)-1 - V;
    constexpr int degree_ret = ffd::set_theory::CardinalitySet(set);
    NilpotentPolynomialF<degree_ret, Field> ret;
    ret.set = set;
    for(BinaryInt S = set; S != 0; S = ((S-1)&set)){
      ret[S] = P_[S];
    }
    ret[0] = P_[0];
    return ret;
  }
  
  
}//namespace ffd::nilpotent_polynomial

