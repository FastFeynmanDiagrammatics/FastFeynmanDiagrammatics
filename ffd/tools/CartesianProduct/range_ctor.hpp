namespace ffd::user_space{

  template<auto d>

  auto

  CartesianProductRange(std::array<std::array<int, 2>, d> L){

    std::array<std::vector<int>, d> ranges;
    
    using ffd::vector_range::Range;
    
    for( int j: Range(d) ){
      ranges[j] = Range(L[j][0], L[j][1]);
    }

    return CartesianProduct<d, std::vector, int>( ranges );
  }

  

  template<auto d>

  auto

  CartesianProductRange(std::array<int, d> L){

    using ffd::vector_range::Range;
    
    std::array<std::array<int, 2>, d> LL;

    for( int j: Range(d) ){
      LL[j][0] = 0;
      LL[j][1] = L[j];
    }
    
    return CartesianProductRange( LL );
  }


  
  template<auto d>

  auto

  CartesianProductRange(int L_low, int L_high){
    
    std::array<std::array<int, 2>, d> LL;
    for( int j: ffd::vector_range::Range(d) ){
      LL[j] = std::array<int, 2>{L_low, L_high};
    }

    return CartesianProductRange( LL );
  }


  
  template<auto d>

  auto

  CartesianProductRange(int L){
    return CartesianProductRange<d>(0, L);
  }




}//namespace
