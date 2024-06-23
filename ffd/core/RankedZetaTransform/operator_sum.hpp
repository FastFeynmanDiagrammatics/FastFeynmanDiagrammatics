namespace ffd::ranked_zeta{

	template<typename Field, SetT SetType>
  RankedZeta<Field, SetType>&
  RankedZeta<Field, SetType>::
  operator+=(const RankedZeta<Field, SetType>& Z) noexcept{
    if(!Z.NotZero){
      return *this;
    }
    std::size_t Z_size = Z.size();
    assert(this->size() == Z_size);
    assert(this->rank_size == Z.rank_size);

    for(size_t S = 0; S < Z_size; ++S)
    {
      (*this)[S] += Z[S];
    }
    return *this;
  }

  template<typename Field, SetT SetType>
  RankedZeta<Field, SetType>&
  RankedZeta<Field, SetType>::
  operator-=(const RankedZeta<Field, SetType>& Z) noexcept{
    return *this += -Z;
  }

  template<typename Field, SetT SetType>
  RankedZeta<Field, SetType>
  RankedZeta<Field, SetType>::operator+(const RankedZeta<Field, SetType>& Z) const noexcept{
    RankedZeta<Field, SetType> ret = *this;
    return ret += Z;
  }

  template<typename Field, SetT SetType>
  RankedZeta<Field, SetType>
  RankedZeta<Field, SetType>::operator-(const RankedZeta<Field, SetType>& add_) const noexcept{
    RankedZeta<Field, SetType> ret = -add_;
    return ret += *this;
  }
}//namespace
