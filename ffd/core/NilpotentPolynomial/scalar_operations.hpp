namespace ffd::nilpotent_polynomial{

  template<typename Field, SetT SetType>

  NilpotentPolynomial<Field, SetType>&

  NilpotentPolynomial<Field, SetType>::operator+=(const Field number_) noexcept{
    (*this)[0] += number_;
    return *this;
  }


  template<typename Field, SetT SetType>
  NilpotentPolynomial<Field, SetType>& NilpotentPolynomial<Field, SetType>::operator-=(const Field number_) noexcept{
    return *this += -number_;
  }


  template<typename Field, SetT SetType>
  NilpotentPolynomial<Field, SetType> NilpotentPolynomial<Field, SetType>::operator+(const Field number_) const noexcept{
    NilpotentPolynomial<Field, SetType> ret = *this;
    return ret += number_;
  }


  template<typename Field, SetT SetType>
  NilpotentPolynomial<Field, SetType> NilpotentPolynomial<Field, SetType>::operator-(const Field number_) const noexcept{
    NilpotentPolynomial<Field, SetType> ret = *this;
    return ret -= number_;
  }


  template<typename Field, SetT SetType>
  NilpotentPolynomial<Field, SetType>& NilpotentPolynomial<Field, SetType>::operator*=(const Field number_) noexcept{
    if(this->NotZero){
      for(BinaryInt S=0; S < this->size(); ++S){
	this->coef[S] *= number_;
      }
    }
    return *this;
  }


  template<typename Field, SetT SetType>
  NilpotentPolynomial<Field, SetType> NilpotentPolynomial<Field, SetType>::operator*(const Field number_) const noexcept{
    NilpotentPolynomial<Field, SetType> ret = *this;
    return ret *= number_;
  }


  template<typename Field, SetT SetType>
  NilpotentPolynomial<Field, SetType>& NilpotentPolynomial<Field, SetType>::operator/=(const Field number_) noexcept{
    *this *= 1./number_;
    return *this;
  }


  template<typename Field, SetT SetType>
  NilpotentPolynomial<Field, SetType> NilpotentPolynomial<Field, SetType>::operator/(const Field number_) const noexcept{
    NilpotentPolynomial<Field, SetType> ret = *this;
    return ret /= number_;
  }


  //with this we can write -P[[\xi]]
  template<typename Field, SetT SetType>
  NilpotentPolynomial<Field, SetType>
  NilpotentPolynomial<Field, SetType>::operator-() const noexcept{
    auto ret = *this;
    if( ret.NotZero ){
      for( std::size_t j = 0ul; j < ret.size(); ++j){
	ret.coef[j] = -ret.coef[j];
      }
    }
    return ret;
  }


  //with this we can write +P[[\xi]]
  template<typename Field, SetT SetType>
  NilpotentPolynomial<Field, SetType>
  NilpotentPolynomial<Field, SetType>::operator+() const noexcept{
    return *this;
  }


  template<typename Field, SetT SetType>
  NilpotentPolynomial<Field, SetType> operator+(Field number_, NilpotentPolynomial<Field, SetType> const& poly_) noexcept{
    return poly_+number_;
  }


  template<typename Field, SetT SetType>
  NilpotentPolynomial<Field, SetType> operator-(Field number_, NilpotentPolynomial<Field, SetType> const& poly_) noexcept{
    return -poly_ + number_;
  }


  template<typename Field, SetT SetType>
  NilpotentPolynomial<Field, SetType> operator*(Field number_, NilpotentPolynomial<Field, SetType> const& poly_) noexcept{
    return poly_ * number_;
  }

  template<typename Field, SetT SetType>
  NilpotentPolynomial<Field, SetType> operator/(Field number_, NilpotentPolynomial<Field, SetType> const& poly_) noexcept{
    NilpotentPolynomial<Field, SetType> ret(__builtin_popcount(poly_.size()-1)); //this is necessary because we have no clue
    ret[0] = number_;
    return ret /= poly_;
  }


}//namespace ffd::nilpotent_polynomial
