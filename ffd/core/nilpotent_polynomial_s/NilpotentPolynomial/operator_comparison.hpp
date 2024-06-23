namespace ffd::nilpotent_polynomial{

	template<typename Field, SetT SetType>
  bool
  NilpotentPolynomial<Field, SetType>::
  operator<(const NilpotentPolynomial<Field, SetType>& P) const noexcept{
    return ((*this)[0] < P[0]);
  }

  template<typename Field, SetT SetType>
  bool
  NilpotentPolynomial<Field, SetType>::
  operator>(const NilpotentPolynomial<Field, SetType>& P) const noexcept{
    return ((*this)[0] > P[0]);
  }

  template<typename Field, SetT SetType>
  bool
  NilpotentPolynomial<Field, SetType>::
  operator<=(const NilpotentPolynomial<Field, SetType>& P) const noexcept{
    return ((*this)[0] <= P[0]);
  }

  template<typename Field, SetT SetType>
  bool
  NilpotentPolynomial<Field, SetType>::
  operator>=(const NilpotentPolynomial<Field, SetType>& P) const noexcept{
    return ((*this)[0] >= P[0]);
  }

  template<typename Field, SetT SetType>
  bool
  NilpotentPolynomial<Field, SetType>::
  operator==(const NilpotentPolynomial<Field, SetType>& P) const noexcept{
    std::size_t const P_size = P.size();
    if (this->size() != P_size)
    {
      return false;
    }
    for(BinaryInt S = 0; S < P_size; ++S){
      if ((*this)[S] != P[S]){
        return false;
      }
    }
    return true;
  }

  template<typename Field, SetT SetType>
  bool
  NilpotentPolynomial<Field, SetType>::
  operator<(const Field F) const noexcept{
    return ((*this)[0] < F);
  }

  template<typename Field, SetT SetType>
  bool
  NilpotentPolynomial<Field, SetType>::
  operator>(const Field F) const noexcept{
    return ((*this)[0] > F);
  }

  template<typename Field, SetT SetType>
  bool
  NilpotentPolynomial<Field, SetType>::
  operator<=(const Field F) const noexcept{
    return ((*this)[0] <= F);
  }

  template<typename Field, SetT SetType>
  bool
  NilpotentPolynomial<Field, SetType>::
  operator>=(const Field F) const noexcept{
    return ((*this)[0] >= F);
  }

  template<typename Field, SetT SetType>
  bool
  operator>=(Field F, const NilpotentPolynomial<Field, SetType> P) noexcept{
	return (F >= P[0]);
  }

  template<typename Field, SetT SetType>
  bool
  operator>(Field F, const NilpotentPolynomial<Field, SetType> P) noexcept{
  return (F > P[0]);
  }

  template<typename Field, SetT SetType>
  bool
  operator<=(Field F, const NilpotentPolynomial<Field, SetType> P) noexcept{
  return (F <= P[0]);
  }

  template<typename Field, SetT SetType>
  bool
  operator<(Field F, const NilpotentPolynomial<Field, SetType> P) noexcept{
  return (F < P[0]);
  }

}//namespace
