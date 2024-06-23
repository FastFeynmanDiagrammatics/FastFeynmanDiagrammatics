namespace ffd::ranked_zeta{

  template<typename Field, SetT SetType>
  [[nodiscard]] RankedZeta<Field, SetType>
  RankedZeta<Field, SetType>::
  extend_to_set(BinaryInt set, int order) const noexcept{

    if(!this->NotZero){
      return RankedZeta<Field, SetType>();
    }

    RankedZeta<Field, SetType> ret(order);

    BinaryInt S = set, S_i = cardinality-1;
    for(; S != 0 && S_i != 0; S = ((S-1)&set), --S_i){
      for (int k=0; k<=(int)rank_size; ++k){
	ret[S+k*ret.cardinality] = (*this)[S_i+k*cardinality];
      }
    }
    ret[0] = (*this)[0];
    return ret;
  }

}//namespace
