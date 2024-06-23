namespace ffd::math_tools{

  template<typename vector_summable_t>
  auto HalfSumDiff(vector_summable_t x){
    return std::array<typename std::decay<decltype(x[0])>::type, 2>
      {.5*(x[0]+x[1]), .5*(x[0]-x[1])};
  }


  template<typename summable_t0, typename summable_t1>
  auto HalfSumDiff(summable_t0 add0, summable_t1 add1){
    std::array<typename std::decay<decltype(add0+add1)>::type, 2> x;
    x[0] = add0; x[1] = add1;
    return HalfSumDiff(x);
  }


}//namespace
