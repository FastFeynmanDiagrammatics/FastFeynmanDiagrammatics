
namespace ffd::nilpotent_polynomial{

  template<typename Field, SetT SetType>
  [[nodiscard]] NilpotentPolynomial<Field, SetType>
  NilpotentPolynomial<Field, SetType>::
  shift(BinaryInt set, int order) const noexcept{
    if( !this->NotZero ){
      return NilpotentPolynomial<Field, SetType>();
    }
    NilpotentPolynomial<Field, SetType> ret(order);


    for( BinaryInt S=0; S < this->coef.size(); ++S){
      if( (S&set) == 0 ){
	if( S+set < ret.coef.size() ){
	  ret[S+set] = (*this)[S];
	}
      }
    }

    return ret;
  }
}//namespace
