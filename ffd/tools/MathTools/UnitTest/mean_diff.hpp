namespace ffd::math_tools::unit_test{

  void mean_diff(){

    [[maybe_unused]] auto [m, d] = HalfSumDiff(std::array{.2, .3});
    [[maybe_unused]] auto [m0, d0] = HalfSumDiff(2, 3.);
    
    // std::cerr<<m<<" "<<d<<std::endl;
    // std::cerr<<m0<<" "<<d0<<std::endl;
    
  }
  

}//namespace
