namespace ffd::large_array::unit_test{

  void test(){
    std::size_t constexpr n = 3;
    l_array<Real, n> v(0);


    std::cerr<<v[0]<<" "<<v[1]<<" "<<v[2]<<" ";
    // std::cerr<<v[3]<<std::endl;

  }

}//namespace
