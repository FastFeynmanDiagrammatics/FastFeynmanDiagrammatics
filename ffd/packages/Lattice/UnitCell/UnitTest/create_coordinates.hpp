

namespace ffd::lattice::unit_test{

  void create_coordinates_unit_cell(){
    using ffd::size, ffd::get, std::size;

    SegmentUnitCell<2> S;
    auto X = CreateCoordinates(S);
    static_assert( size(X) == 1 );
    static_assert( std::is_same<decltype(X),
		   ffd::periodic_coordinate::PeriodicCoordinates<2, int>>::value );

    assert( get<0>(X)() == 0 );
    get<0>(X).Variable = 3;
    assert( get<0>(X)() == 1 );


    EquilateralTriangleUnitCell<3> T;
    auto Y = CreateCoordinates(T);
    static_assert( size(Y) == 1 );
    static_assert( std::is_same<decltype(Y),
		   ffd::periodic_coordinate::PeriodicCoordinates<3, int>>::value );

    assert( get<0>(Y)() == 0 );
    get<0>(Y).Variable = 3;
    assert( get<0>(Y)() == 0 );
    get<0>(Y).Variable = 4;
    assert( get<0>(Y)() == 1 );

    assert( size(get<0>(Y).RealSpaceVectors) == 3);
    std::array<std::array<Real, 3>, 3> real_space_vectors;
    real_space_vectors[0] = {0.};
    real_space_vectors[1] = {1., 0., 0.};
    real_space_vectors[2] = {.5, .5*std::sqrt(3.), 0.};
    for(auto j: {0, 1, 2}){
      for(auto k: {0, 1, 2}){
	assert( std::abs(real_space_vectors[j][k] - get<0>(Y).RealSpaceVectors[j][k]) <
		std::numeric_limits<Real>::epsilon() );
      }
    }
    
  }

}//namespace
