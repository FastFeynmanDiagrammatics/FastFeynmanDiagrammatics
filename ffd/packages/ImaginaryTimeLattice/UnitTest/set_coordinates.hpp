namespace ffd::imaginary_time_lattice::unit_test{

  void set_coordinates(){

    Real beta = 2.1232;
    

    ffd::imaginary_time::PeriodicImaginaryTime one_over_T(beta);
    ffd::lattice::HoneycombLattice H({10, 12});


    Real tau = .234678324798;
    std::array<int, 3> x_set;
    x_set[0] = 4;
    x_set[1] = 2;
    x_set[2] = 3;
    auto X = SetCoordinates(one_over_T, H, std::make_pair(tau, x_set));
    

    assert(( std::abs( ffd::get<0>(X)()-tau) < 2*std::numeric_limits<Real>::epsilon() ));
    assert(( ffd::get<1>(X)() == 4 ));
    assert(( ffd::get<2>(X)() == 2 ));
    assert(( ffd::get<3>(X)() == 1 )); 
    

  }

}//namespace
