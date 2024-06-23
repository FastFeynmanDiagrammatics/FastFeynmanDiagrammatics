

namespace ffd::lattice::unit_test{

  void triangular_lattice(){
    TriangularLattice T({4, 4});

    for(int j=0; j < 2; ++j){
      assert( T.LowerUpperBounds[j][0] == -1);
      assert( T.LowerUpperBounds[j][1] ==  3);
    }
    std::array<std::array<Real, 2>, 2> bravais_vectors;
    bravais_vectors[0] = {1., 0.};
    bravais_vectors[1] = {.5, .5*std::sqrt(3.)};
    for(int j=0; j<2; ++j){
      for(int k=0; k<2; ++k){
	assert( std::abs( bravais_vectors[j][k] - T.BravaisVectors[j][k]) <
		10*std::numeric_limits<Real>::epsilon());
      }
    }


    TriangularLattice T_half_pi({3, 3}, 1., .5);

    std::array<std::array<Real, 2>, 2> bravais_vectors_half_pi;
    bravais_vectors_half_pi[0] = {0., 1.};
    bravais_vectors_half_pi[1] = {-.5*std::sqrt(3.), .5};
    for(int j=0; j<2; ++j){
      for(int k=0; k<2; ++k){
	assert( std::abs( bravais_vectors_half_pi[j][k] - T_half_pi.BravaisVectors[j][k]) <
		10*std::numeric_limits<Real>::epsilon());
      }
    }


    TriangularLattice T_sixth_pi({3, 3}, 1., 1./6.);

    std::array<std::array<Real, 2>, 2> bravais_vectors_sixth_pi;
    bravais_vectors_sixth_pi[0] = {.5*sqrt(3.), .5};
    bravais_vectors_sixth_pi[1] = {0., 1.};
    for(int j=0; j<2; ++j){
      for(int k=0; k<2; ++k){
	assert( std::abs( bravais_vectors_sixth_pi[j][k] - T_sixth_pi.BravaisVectors[j][k]) <
		10*std::numeric_limits<Real>::epsilon());
      }
    }

  }

}//namespace
