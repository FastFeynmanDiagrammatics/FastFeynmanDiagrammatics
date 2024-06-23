

namespace ffd::lattice::unit_test{
  using namespace ffd::user_space;
  void create_coordinates(){
    using ffd::size, ffd::get;
    TriangularLattice T({3, 3});

    auto X = CreateCoordinates(T);
    static_assert( size(X) == 2);
    static_assert( std::is_same<decltype(X),
		   ffd::periodic_coordinate::PeriodicCoordinates<2, int, int>>::value );

    std::array<std::array<int, 2>, 2> lower_upper_bounds;
    for(auto j: {0, 1}){
      lower_upper_bounds[j] = {-1, 2};
    }
    for(auto j: {0, 1}){
      for(auto k: {0, 1}){
	assert( std::abs(lower_upper_bounds[j][k] - T.LowerUpperBounds[j][k]) <
		10*std::numeric_limits<Real>::epsilon() );
      }
    }

    std::array<std::vector<std::array<Real, 2>>, 2> real_space_vectors;
    real_space_vectors[0].push_back({1., 0});
    real_space_vectors[1].push_back({.5, .5*sqrt(3.)});
    for(auto k: {0, 1}){
      assert( std::abs(real_space_vectors[0][0][k] - get<0>(X).RealSpaceVectors[0][k])
	      < 10*std::numeric_limits<Real>::epsilon() );
    }
    for(auto k: {0, 1}){
      assert( std::abs(real_space_vectors[1][0][k] - get<1>(X).RealSpaceVectors[0][k])
	      < 10*std::numeric_limits<Real>::epsilon() );
    }

  }

}//namespace
