

namespace ffd::nilpotent_polynomial{
  
  //for the moment assume nilpotency == 1 ( != 1 not correctly implemented)
  //TODO implement nilpotency != 1
  template <int degree, typename Field = Real, int nilpotency = 1>
  struct NilpotentPolynomialR{
    //the numeric base we use for (multi)set-theory computations.
    //for standard set theory, we use base 2
    static const int NumericBase = nilpotency + 1;

    //the number of coefficients of the polynomials
    static const BinaryInt NumberCoef = ffd::core_math::Pow1ToInt2(NumericBase, degree);

    static const BinaryInt NumberSubsets = ffd::core_math::Pow1ToInt2(NumericBase+1, degree);

    //this is used just for printing
    static const char variable_names[26];
  
    //Sets and Subsets are used for multiplications and division
    //(see documentation)
    static BinaryInt Sets[NumberSubsets-NumberCoef+(degree==0)];
    static BinaryInt Subsets[NumberSubsets-NumberCoef+(degree==0)];
  
    //this is equal to 1 if Sets and Subsets have never been computed
    static bool ClassNotInitialized;

    //the subset of \{\xi_1, ..., \xi_n\} in which we perform operations
    BinaryInt set = NumberCoef-1;

    //the array of coefficients
    Field coef[NumberCoef] = {0.};

    //this bool answer the question whether this polynomial is the zero polynomial
    bool ThisIsNotTheZeroPolynomial = true;
  
    NilpotentPolynomialR() {}
    
    NilpotentPolynomialR(const char* WhichPolynomial){
      //std::cerr<<"STRN";
      for(BinaryInt V=0; V<NumberCoef; ++V){
	coef[V] = 0;
      }
      set = 0;
      if(WhichPolynomial == "Zero"){
	ThisIsNotTheZeroPolynomial = false;
      }else if(WhichPolynomial == "One"){
	coef[0] = 1.;
      }else{
	throw std::invalid_argument("NilpotentPolynomialR::NilpotentPolynomialR(std::string): argument unknown");
      }
    }

    
    NilpotentPolynomialR(Field x){
      //std::cerr<<"FILD";
      for(BinaryInt V=1; V<NumberCoef; ++V){
	coef[V] = 0.;
      }
      coef[0] = x;
      set = 0;
    }

    
    NilpotentPolynomialR& operator=(const NilpotentPolynomialR& source){
      for(BinaryInt V=0; V<NumberCoef; ++V){
    	coef[V] = source.coef[V];
      }
      set = source.set;
      ThisIsNotTheZeroPolynomial = source.ThisIsNotTheZeroPolynomial;
      //std::cerr<<"EQUL";
      return *this;
    }
    
    
    NilpotentPolynomialR(const NilpotentPolynomialR& source){
      //std::cerr<<"COPY";
      *this = source;
    }
    
  
    //returns the maximal value of coef
    Real MaxAbsCoef(){
      Real abs_err_max = 0;
      for(BinaryInt V=0; V<NumberCoef; ++V){
	if(abs_err_max < std::abs(coef[V]))
	  abs_err_max = std::abs(coef[V]);
      }
      return abs_err_max;
    }
  
  
  
    //Initialize the Sets and Subsets static arrays shared by all objects of the class.
    //For this reason, it should be called only once. If the class was already initialized
    //in the past, the bool ClassNotInitialized is false.
    void InitializeClass(){
      ClassNotInitialized = false;
      BinaryInt array_index = 0;
      for(BinaryInt V=0; V<NumberCoef; ++V){
	for(BinaryInt S=0; S<V; ++S){
	  if(ffd::set_theory::IsSet1SubsetOfSet2(S, V, nilpotency)){
	    Sets[array_index] = V;
	    Subsets[array_index] = S;
	    ++array_index;
	  }
	}
      }
    }



    //it converts a set S belonging to the internal convention
    //(which means that S \in {0, 1, 2^degree-1})
    //in a set belonging to the external convention,
    //which is determined by "set"
    inline BinaryInt ConvertInternalSetToExternal(BinaryInt S) const{
      BinaryInt set_copy = set, V = 0;
      for(BinaryInt pow_NumericBase = 1; set_copy>0; pow_NumericBase *= NumericBase){
	if(set_copy%NumericBase > 0){
	  V += (S%NumericBase)*pow_NumericBase;
	  S /= NumericBase;
	}
	set_copy /= NumericBase;
      }
      return V;
    }

  
  
    //it converts a set belonging to the external convention
    //to a set belonging to the internal convention
    //for coefficients
    BinaryInt ConvertExternalSetToInternal(BinaryInt S){
      BinaryInt set_copy = set, V = 0;
      for(BinaryInt pow_NumericBase = 1; S>0; S /= NumericBase){
	V += (S%NumericBase)*pow_NumericBase;
	if(set_copy%NumericBase > 0){
	  pow_NumericBase *= NumericBase;
	}
	set_copy/= NumericBase;
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

    NilpotentPolynomialR& operator+=(Field number){
      ThisIsNotTheZeroPolynomial = true;
      coef[0] += number;
      return *this;
    }

    NilpotentPolynomialR& operator-=(Field number){
      return *this += -number;
    }
    
    NilpotentPolynomialR operator+(Field number) const{
      NilpotentPolynomialR ret = *this;
      return ret += number;
    }

    
    NilpotentPolynomialR operator-(Field number) const{
      NilpotentPolynomialR ret = *this;
      return ret -= number;
    }

    
    NilpotentPolynomialR& operator*=(Field number){
      if(ThisIsNotTheZeroPolynomial){
	for(BinaryInt S=0; S<NumberCoef; S++){
	  coef[S] *= number;
	}
      }
      return *this;
    }

    
    NilpotentPolynomialR operator*(Field number) const{
      NilpotentPolynomialR ret = *this;
      return ret *= number;
    }

    
    //with this we can write -P[[\xi]]
    NilpotentPolynomialR operator-() const{
      NilpotentPolynomialR ret = *this;
      return ret *= -1.;
    }

    NilpotentPolynomialR& operator/=(Field number){
      *this *= 1./number;
      return *this;
    }

    NilpotentPolynomialR operator/(Field number) const{
      NilpotentPolynomialR ret = *this;
      return ret /= number;
    }
    
    
    //The vector product multiply each \xi_j by some vector component v_j  
    NilpotentPolynomialR operator*(std::vector<Field> weight) const{
      NilpotentPolynomialR ret = *this;
      for(BinaryInt S=0; S<NumberCoef; S++){
	for(int j=0; j<degree; ++j){
	  BinaryInt pow_NumericBase = ffd::core_math::Pow1ToInt2(NumericBase, j);
	  for(int l=1; l<NumericBase; ++l){
	    if(ffd::set_theory::IsSet1SubsetOfSet2(pow_NumericBase*l, ConvertInternalSetToExternal(S), nilpotency)){
	      ret[S] *= weight[j];
	    }
	  }
	}
      }
    }

    
    NilpotentPolynomialR& operator*=(std::vector<Field> weight){
      return *this = *this * weight;
    }

    
    NilpotentPolynomialR& operator+=(const NilpotentPolynomialR& add){
      if(add.ThisIsNotTheZeroPolynomial){
	ThisIsNotTheZeroPolynomial = true;
      }
      coef[0] += add[0];
      assert(set==0 || add.set==0 || set == add.set);
      if(add.set != 0){
	for(BinaryInt S=1; S<NumberCoef; ++S){
	  coef[S] += add[S];
	}
      }
      if(set==0){
	set = add.set;
      }
      return *this;
    }


    NilpotentPolynomialR& operator-=(const NilpotentPolynomialR& add){
      return *this += -add;
    }

    
    //addition between polynomials
    NilpotentPolynomialR operator+(const NilpotentPolynomialR& A1) const{
      NilpotentPolynomialR ret = *this;
      return ret += A1;
    }

    NilpotentPolynomialR operator-(const NilpotentPolynomialR& A1) const{
      NilpotentPolynomialR ret = -A1;
      return ret += *this;
    }

    //multiply *this by a polynomial 
    NilpotentPolynomialR operator*(const NilpotentPolynomialR& fac) const{
      if(!ThisIsNotTheZeroPolynomial || !fac.ThisIsNotTheZeroPolynomial){
	return "Zero";
      }
      assert(set==0 || fac.set==0 || set == fac.set);
      NilpotentPolynomialR ret;
      if(ClassNotInitialized){
	ret.InitializeClass();
      }
      if(set == 0){
	if(fac.set == 0){
	  ret[0] = coef[0]*fac[0];
	}else{
	  ret = fac*coef[0];
	}
      }else{
	ret = *this * fac[0];
	ret.set = set;
	if(fac.set != 0){
	  for(BinaryInt V=0; V<NumberSubsets-NumberCoef; ++V){
	    ret[Sets[V]] += coef[Subsets[V]]*fac[Sets[V]-Subsets[V]];
	  }
	}
      }
      return ret;
    }
    
    
    //multiplication between polynomials
    NilpotentPolynomialR& operator*=(const NilpotentPolynomialR& fac){
      return *this = *this * fac;
    }

    
    //divide *this by a polynomial
    NilpotentPolynomialR operator/(NilpotentPolynomialR den) const{
      if(!ThisIsNotTheZeroPolynomial){
	return *this;
      }
      assert(set==0 || den.set==0 || set == den.set);
      NilpotentPolynomialR ret;
      if(ClassNotInitialized){
	ret.InitializeClass();
      }
      if(den.set == 0){
	if(set == 0){
	  ret[0] = coef[0]/den[0];
	}
	else{
	  ret = *this/den[0];
	}
      }
      else{
	ret.set = den.set;
	Field den_empty_set = den[0];
	NilpotentPolynomialR num = *this;
	num /= den_empty_set;
	ret = num;
	den /= den_empty_set;
	for(BinaryInt V=0; V < NumberSubsets - NumberCoef; ++V){
	  ret[Sets[V]] -= ret[Subsets[V]]*den[Sets[V]-Subsets[V]];
	}
      }
      return ret;
    }


    NilpotentPolynomialR& operator/=(const NilpotentPolynomialR& den){
      return *this = *this / den;
    }
    
    
    //restrict the polynomial over a specified set S (that can be of lower order r)
    //and it stores it in the second argument
    //WARNING: it works only if nilpotency == 1
    template<int r>
    void RestrictTo1StoreIn2(BinaryInt S, NilpotentPolynomialR<r, Field, nilpotency>& dest){
      if(!ThisIsNotTheZeroPolynomial){
	dest.ThisIsNotTheZeroPolynomial = false;
      }
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



  template <int degree, typename Field, int nilpotency>
  const char NilpotentPolynomialR<degree, Field, nilpotency>::variable_names[] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'};


  template <int degree, typename Field, int nilpotency>
  bool NilpotentPolynomialR<degree, Field, nilpotency>::ClassNotInitialized = true;
  

  template<int degree, typename Field, int nilpotency>
  BinaryInt NilpotentPolynomialR<degree, Field, nilpotency>::Sets[] = {0};
  

  template<int degree, typename Field, int nilpotency>
  BinaryInt NilpotentPolynomialR<degree, Field, nilpotency>::Subsets[] = {0};
  
}//namespace ffd::nilpotent_polynomial


