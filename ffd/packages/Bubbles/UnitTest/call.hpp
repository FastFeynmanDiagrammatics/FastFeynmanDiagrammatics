namespace ffd::bubbles::unit_test{

  void UnitTest(){

    assert(( gamma0_atom(2.221412, 1.123412, 1.423879) ));

    assert(( tadpoles_atom(2.221412, 1.123412, 1.423879) ));

    assert(( tadpoles_atom(1, -.311282, .1) ));



    int L = 10;
    Real beta = 10;
    Real mu0 = -.3;
    gamma0_square_lattice<true, true>(beta, -.3, 2, L);
    
  }

}//namespace
