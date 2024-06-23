namespace ffd::nilpotent_polynomial{

  template<typename Field, SetT SetType>
  [[nodiscard]] NilpotentPolynomial<Field, SetType>
  NilpotentPolynomial<Field, SetType>::
  extend_to_set(BinaryInt set, int order) const noexcept{

    if(!this->NotZero){
      return NilpotentPolynomial<Field, SetType>();
    }
    int const order_now = this->linear_size();

    NilpotentPolynomial<Field, SetType> ret(order);

    BinaryInt S = set, S_i = (1<<order_now)-1;
    for(; S != 0 && S_i != 0; S = ((S-1)&set), --S_i){
      ret[S] = (*this)[S_i];
    }
    ret[0] = (*this)[0];
    return ret;
  }

}//namespace
