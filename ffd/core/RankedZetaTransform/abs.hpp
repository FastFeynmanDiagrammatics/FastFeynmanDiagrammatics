namespace ffd::ranked_zeta{

	template<typename Field, SetT SetType>
    RankedZeta<Field, SetType>
    abs(const RankedZeta<Field, SetType>& Z){
	  int order = Z.order();
      RankedZeta<Field, SetType> Z_abs(order);
      for( std::size_t j = 0ul; j < Z_abs.size(); ++j){
        Z_abs[j] = std::abs(Z[j]);
      }
      return Z_abs;
    }

}//namespace
