namespace ffd::ranked_zeta{

  template<typename Field, SetT SetType>
  bool
  RankedZeta<Field, SetType>::
  operator<(const RankedZeta<Field, SetType>& P) const noexcept{
    return ((*this)[0] < P[0]);
  }

  template<typename Field, SetT SetType>
  bool
  RankedZeta<Field, SetType>::
  operator>(const RankedZeta<Field, SetType>& P) const noexcept{
    return ((*this)[0] > P[0]);
  }

  template<typename Field, SetT SetType>
  bool
  RankedZeta<Field, SetType>::
  operator<=(const RankedZeta<Field, SetType>& P) const noexcept{
    return ((*this)[0] <= P[0]);
  }

  template<typename Field, SetT SetType>
  bool
  RankedZeta<Field, SetType>::
  operator>=(const RankedZeta<Field, SetType>& P) const noexcept{
    return ((*this)[0] >= P[0]);
  }

  template<typename Field, SetT SetType>
  bool
  RankedZeta<Field, SetType>::
  operator==(const RankedZeta<Field, SetType>& Z) const noexcept{
    std::size_t const Z_size = Z.size();
    if (this->size() != Z_size)
    {
      return false;
    }
    for(BinaryInt S = 0; S < Z_size; ++S){
      if ((*this)[S] != Z[S]){
        return false;
      }
    }
    return true;
  }

  template<typename Field, SetT SetType>
  bool
  RankedZeta<Field, SetType>::
  operator<(const Field F) const noexcept{
    return ((*this)[0] < F);
  }

  template<typename Field, SetT SetType>
  bool
  RankedZeta<Field, SetType>::
  operator>(const Field F) const noexcept{
    return ((*this)[0] > F);
  }

  template<typename Field, SetT SetType>
  bool
  RankedZeta<Field, SetType>::
  operator<=(const Field F) const noexcept{
    return ((*this)[0] <= F);
  }

  template<typename Field, SetT SetType>
  bool
  RankedZeta<Field, SetType>::
  operator>=(const Field F) const noexcept{
    return ((*this)[0] >= F);
  }

  template<typename Field, SetT SetType>
  bool
  RankedZeta<Field, SetType>::
  operator==(const Field F) const noexcept{
    return ((*this)[0] == F);
  }

  template<typename Field, SetT SetType>
  bool
  operator>=(Field F, const RankedZeta<Field, SetType> Z) noexcept{
  return (F >= Z[0]);
  }

  template<typename Field, SetT SetType>
  bool
  operator>(Field F, const RankedZeta<Field, SetType> Z) noexcept{
  return (F > Z[0]);
  }

  template<typename Field, SetT SetType>
  bool
  operator<=(Field F, const RankedZeta<Field, SetType> Z) noexcept{
  return (F <= Z[0]);
  }

  template<typename Field, SetT SetType>
  bool
  operator<(Field F, const RankedZeta<Field, SetType> Z) noexcept{
  return (F < Z[0]);
  }

}//namespace
