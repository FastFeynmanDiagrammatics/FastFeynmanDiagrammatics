

namespace ffd::user_space::imaginary_time_lattice_propagator::unit_test{

  void test_get_bravais_vectors(){
    int L = 5;
    ffd::lattice::HoneycombLattice H{{L, L}};

    std::vector<Real> bravais_vectors(4, 0.);
    auto x = CreateCoordinates(H);
    get_bravais_vectors<2, 0>(x, bravais_vectors);

  }


}//namespace
