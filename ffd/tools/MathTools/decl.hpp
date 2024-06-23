namespace ffd::math_tools{

  template<typename IntType>
  constexpr int log2_int(IntType x){

    static_assert( std::is_integral_v<IntType> );
    
    unsigned long log2 = 0;
    for(;; ++log2){
      if( (1ul<<log2) > (unsigned long)(x)){
	--log2;
	break;
      }
    }

    return log2;
  }

}//namespace


