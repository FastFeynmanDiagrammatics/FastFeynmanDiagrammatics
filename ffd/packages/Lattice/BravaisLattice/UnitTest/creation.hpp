

namespace ffd::lattice::unit_test{

  void creation(){
    
    BravaisLattice<1> L({1});

    std::array<std::array<int, 2>, 1> vec;
    vec[0] = {0, 1};
    assert(L.LowerUpperBounds == vec);


    BravaisLattice<1> L2({2});

    vec[0] = {0, 2};
    assert(L2.LowerUpperBounds == vec);


    BravaisLattice<1> L3({3});

    vec[0] = {-1, 2};
    assert(L3.LowerUpperBounds == vec);


    BravaisLattice<1> L4({4});

    vec[0] = {-1, 3};
    assert(L4.LowerUpperBounds == vec);

  }

}//namespace
