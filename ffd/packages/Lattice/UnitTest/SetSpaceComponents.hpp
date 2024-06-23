namespace ffd::lattice::unit_test{

  void set_space_components(){

    int const Lx = 3, Ly = 5;

    
    HypercubicLattice<2> H({Lx, Ly});


    
    std::array<int, 3> x{1, -2, 2};
      

    auto r = SetSpaceCoordinates(H, x);
    


    assert(( ffd::get<0>(r)() == 1 ));

    assert(( ffd::get<1>(r)() == -2 ));

    assert(( component<0>(r) == 1 ));



    
    KagomeLattice K({Lx, Ly});
    auto r_k = SetSpaceCoordinates(K, x);

    assert(( ffd::get<0>(r_k)() == 1 ));

    assert(( ffd::get<1>(r_k)() == -2 ));

    assert(( ffd::get<2>(r_k)() == 2 ));
    
    
    
  }

}//namespace
