namespace ffd::ranked_zeta{

	template<typename Field, SetT SetType>
    RankedZeta<Field, SetType>
    RankedZeta<Field, SetType>::
    restrict_to_set(BinaryInt set) const noexcept{

      int const order = ffd::set_theory::CardinalitySet(set);
	  BinaryInt const card = this->cardinality;
	  BinaryInt const card_new = (1<<order);

      RankedZeta<Field, SetType> ret(order);

      if( this->size() > 1){
        BinaryInt V_i = card_new-1;
  	    BinaryInt V = set;
        for(; V != 0; V = ((V-1)&set), --V_i){
		   for (int k=0; k<=order; ++k){
  	          ret[V_i+k*card_new] = this->coef[V+k*card];
            }
        }
      }
      for (int k=0; k<=order; ++k){
        ret[k*card_new] = this->coef[k*cardinality];
      }
      ret.NotZero = this->NotZero;
      return ret;
    }

}//namespace
