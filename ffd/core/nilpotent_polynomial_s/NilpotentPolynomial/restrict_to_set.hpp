namespace ffd::nilpotent_polynomial{

  template<typename Field, SetT SetType>
  NilpotentPolynomial<Field, SetType>
  NilpotentPolynomial<Field, SetType>::
  restrict_to_set(BinaryInt set) const noexcept{
    int const order = ffd::set_theory::CardinalitySet(set);
    NilpotentPolynomial<Field, SetType> ret(order);
    if(this->size() > 1){
      BinaryInt V_i = (1<<order)-1,
	  V = set;
      for(; V != 0; V = ((V-1)&set), --V_i){
	    ret[V_i] = (*this)[V];
      }
    }
    ret[0] = (*this)[0];
    ret.NotZero = this->NotZero;
    return ret;
  }



}//namespace ffd::nilpotent_polynomial
