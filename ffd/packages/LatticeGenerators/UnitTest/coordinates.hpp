namespace ffd::lattice_g::unit_test{

  void coordinates(){
    using namespace ffd::user_space;
    auto c = LatticeCoordinates_g<2>();
    auto c12 = c(1, 2);
    std::cerr << c12 << '\n';
  }

}//namespace
