namespace ffd::nilpotent_polynomial{

  template<typename Field, SetT SetType>
  NilpotentPolynomial<Field, SetType>&
  NilpotentPolynomial<Field, SetType>::
  operator+=(const NilpotentPolynomial<Field, SetType>& P) noexcept{
    if(!P.NotZero){
      return *this;
    }
    std::size_t const P_size = P.size();
    if(this->size() < P_size){
      this->coef.resize(P_size);
    }
    for(BinaryInt S = 0; S < P_size; ++S){
      (*this)[S] += P[S];
    }
    return *this;
  }


  template<typename Field, SetT SetType>
  NilpotentPolynomial<Field, SetType>&
  NilpotentPolynomial<Field, SetType>::
  operator-=(const NilpotentPolynomial<Field, SetType>& P) noexcept{
    return *this += -P;
  }


  template<typename Field, SetT SetType>
  NilpotentPolynomial<Field, SetType>
  NilpotentPolynomial<Field, SetType>::operator+(const NilpotentPolynomial<Field, SetType>& P) const noexcept{
    NilpotentPolynomial<Field, SetType> ret = *this;
    return ret += P;
  }


  template<typename Field, SetT SetType>
  NilpotentPolynomial<Field, SetType>
  NilpotentPolynomial<Field, SetType>::operator-(const NilpotentPolynomial<Field, SetType>& add_) const noexcept{
    NilpotentPolynomial<Field, SetType> ret = -add_;
    return ret += *this;
  }


}//namespace ffd::nilpotent_polynomial
