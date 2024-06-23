

namespace ffd::lattice::unit_test{
  
  using namespace ffd::user_space;
  
  void create_coordinates_honeycomb(){
    using namespace std;
    using ffd::get, ffd::RealSpace;
    HoneycombLattice H({10, 10});

    auto X = CreateCoordinates(H);
    static_assert( is_same<decltype(X),
		   ffd::periodic_coordinate::PeriodicCoordinates<2, int, int, int>>::value );
    get<2>(X).Variable = 2;
    assert( get<2>(X)() == 0 );
    assert( get<1>(X)() == 0 );
    assert( get<0>(X)() == 0 );

    auto r = RealSpace(X);
    assert( size(r) == 2);
    for(int j=0; j<2; ++j){
      assert( abs(r[j]) < 10*numeric_limits<Real>::epsilon() );
    }

    get<2>(X).Variable = 1;
    r = RealSpace(X);
    for(int j=0; j<2; ++j){
      assert( abs(r[j] - (j==0)) < 10*numeric_limits<Real>::epsilon() );
    }

    get<1>(X).Variable = 1;
    r = RealSpace(X);
    for(int j=0; j<2; ++j){
      assert( abs(r[j] - (j==0) - (j==1)*sqrt(3.)) < 10*numeric_limits<Real>::epsilon() );
    }

    get<0>(X).Variable = 1;
    auto r2 = RealSpace(X);
    for(int j=0; j<2; ++j){
      assert( abs(r2[j] - r[j] - .5*sqrt(3)*(j==1)-1.5*(j==0)) < 10*numeric_limits<Real>::epsilon() );
    }

    
  }

}//namespace
