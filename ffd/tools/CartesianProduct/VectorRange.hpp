namespace ffd::user_space{

  template<auto d,
  	   typename vector_t>

  auto
  VectorRange(vector_t L){
    std::array<std::array<int, 2>, d> L_array;
    for( int j: ffd::vector_range::Range(d) ){
      L_array[j][0] = 0;
      L_array[j][1] = L[j];
    }
    return CartesianProductRange( L_array );
  }

}//namespace
