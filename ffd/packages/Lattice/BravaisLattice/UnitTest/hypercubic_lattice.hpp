

namespace ffd::lattice::unit_test{

  void hypercubic_lattice(){


    HypercubicLattice<1> L1({1});

    std::array<std::array<Real, 1>, 1> vec1;
    vec1[0] = {1};
    assert(L1.BravaisVectors == vec1);


    HypercubicLattice<2> L2({3, 3});

    std::array<std::array<Real, 2>, 2> vec2;
    vec2[0] = {1, 0};
    vec2[1] = {0, 1};
    assert(L2.BravaisVectors == vec2);


    HypercubicLattice<3> L3({3, 3, 3}, 2.);

    for(int j=0; j < 3; ++j){
      assert(L3.LowerUpperBounds[j][0] == -1);
      assert(L3.LowerUpperBounds[j][1] == 2);
    }
    std::array<std::array<Real, 3>, 3> vec3;
    vec3[0] = {2., 0., 0.};
    vec3[1] = {0., 2., 0.};
    vec3[2] = {0., 0., 2.};
    for(int j=0; j < 3; ++j){
      for(int k=0; k < 3; ++k){
	assert( std::abs(L3.BravaisVectors[j][k]-vec3[j][k]) < 10*std::numeric_limits<Real>::epsilon() );
      }
    }

  }

}//namespace
