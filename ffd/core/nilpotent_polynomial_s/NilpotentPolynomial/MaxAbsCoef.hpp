namespace ffd::nilpotent_polynomial{

  //returns the maximal value of coef
  template<typename Field, SetT SetType>
  Real NilpotentPolynomial<Field, SetType>::MaxAbsCoef() noexcept{
    Real abs_err_max = 0;
    for(BinaryInt V=0; V<coef.size(); ++V){
      if(abs_err_max < std::abs(coef[V])){
	abs_err_max = std::abs(coef[V]);
      }
    }
    return abs_err_max;
  }

}//namespace ffd::nilpotent_polynomial
