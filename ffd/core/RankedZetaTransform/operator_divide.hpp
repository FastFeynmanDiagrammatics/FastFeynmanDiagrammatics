namespace ffd::ranked_zeta{

template<typename Field, SetT SetType>
[[nodiscard]] RankedZeta<Field, SetType>
RankedZeta<Field, SetType>::
operator/(RankedZeta<Field, SetType> const & Z) const noexcept{
  if(!this->NotZero || !Z.NotZero){
    return RankedZeta<Field, SetType>();
  }
  else if(Z.size() == 1){
	return *this/Z[0];
  }
  else if(this->size() == 1){
  return Z/this->coef[0];
  }

  assert(this->size() == Z.size());
  assert(rank_size == Z.rank_size);

  RankedZeta<Field, SetType> ret = *this/Z[0];

  auto z = Z/Z[0];

  BinaryInt A_=0, B_=0, C_=0;

  for (BinaryInt k=0; k<=rank_size; ++k){
    for (BinaryInt j=1; j<=k; ++j){
      A_ += cardinality;
      B_ -= cardinality;
      for (BinaryInt set=0; set<cardinality; ++set)
      {
		ret[C_+set] -= z[A_+set] * ret[B_+set];
      }
    }
    A_ -= k*cardinality;
    B_ += (k+1)*cardinality;
    C_ += cardinality;
  }
  return ret;
  }

  template<typename Field, SetT SetType>
  RankedZeta<Field, SetType>&
  RankedZeta<Field, SetType>::operator/=(const RankedZeta<Field, SetType>& Z) noexcept{
    return *this = *this / Z;
  }
}//namespace
