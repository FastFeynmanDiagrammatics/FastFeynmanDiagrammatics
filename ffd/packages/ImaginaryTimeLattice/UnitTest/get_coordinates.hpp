namespace ffd::imaginary_time_lattice::unit_test{


  void get_coordinates(){

    
    Real beta = 2.1232;
    

    ffd::imaginary_time::PeriodicImaginaryTime one_over_T(beta);
    ffd::lattice::HoneycombLattice H({10, 12});


    Real tau = .234678324798;
    std::array<int, 3> x_set;
    x_set[0] = -4;
    x_set[1] = -2;
    x_set[2] = -3;
    auto Y = SetCoordinates(one_over_T, H, std::make_pair(tau, x_set));
    auto X = GetCoordinates<decltype(one_over_T), decltype(H)>(Y);
    

    assert(( std::abs( X.first - tau) < 2*std::numeric_limits<Real>::epsilon() ));
    assert(( X.second[0] == -4 ));
    assert(( X.second[1] == -2 ));
    assert(( X.second[2] == 1 ));
    // std::cerr<<X.first<<std::endl;
    // std::cerr<<X.second[0]<<std::endl;
    // std::cerr<<X.second[1]<<std::endl;
    // std::cerr<<X.second[2]<<std::endl;

  }    



}//namespace
