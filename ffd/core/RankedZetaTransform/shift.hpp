namespace ffd::ranked_zeta{

	template<typename Field, SetT SetType>
    [[nodiscard]] RankedZeta<Field, SetType>
    RankedZeta<Field, SetType>::
    shift(BinaryInt shift_set, int order) const noexcept{
      if( !this->NotZero ){
        return RankedZeta<Field, SetType>();
      }

	  assert(order <= __builtin_popcount(this->cardinality-1));

	  RankedZeta<Field, SetType> ret(order);
		if (SetType == Even){
			ret.resize((1<<order)*(order/2+1));
			ret.rank_size = order/2;
		}

	  BinaryInt pop_sh = __builtin_popcount(shift_set);
	  if (SetType == Even) {
		  assert(pop_sh%2==0);
		  pop_sh/=2;
	  }

	  for (BinaryInt k=0; k<=ret.rank_size-pop_sh; ++k)
	  {
		  for (BinaryInt set=shift_set; set<(1<<order); ++set)
		  {
			  if((set & shift_set) == shift_set){
				   ret[(k+pop_sh)*ret.cardinality+set] += this->coef[k*cardinality+set-shift_set];
			   }
		  }
	  }
	  // std::cerr<<std::endl;
	  return ret;

    }

}//namespace
