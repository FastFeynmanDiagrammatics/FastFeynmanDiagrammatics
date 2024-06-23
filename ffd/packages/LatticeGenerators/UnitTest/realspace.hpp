namespace ffd::lattice_g::unit_test{

  void realspace(){
    std::array<std::array<Real, 1>, 1> bravais;
    bravais[0] = {1.};
    auto Realspace = Realspace_g<1>({10}, bravais);
    //    std::cerr << Realspace({2})[0] << '\n';
  }

}//namespace
