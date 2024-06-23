namespace ffd::math_tools{

  template<typename IntType1, typename IntType2 = IntType1>

  auto
  DivMod(IntType1 x, IntType2 y){
    std::array<IntType1, 2> ret;

    
    ret[0] = x/y;
    ret[1] = x%y;

    
    return ret;
  }


  template<typename IntType1, typename IntType2 = IntType1>

  auto
  ModDiv(IntType1 x, IntType2 y){
    std::array<IntType1, 2> ret;

    
    ret[1] = x/y;
    ret[0] = x%y;

    
    return ret;
  }

}//namespace

