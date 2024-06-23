

namespace ffd::lattice::unit_test{
  using namespace ffd::user_space;
  
  template<unsigned long size>
  void check_equality_vectors(std::array<Real, size> x1, std::array<Real, size> x2){
    using std::size;
    for(std::size_t j=0; j < size(x1); ++j){
      assert( std::abs(x1[j]-x2[j]) < 10*std::numeric_limits<Real>::epsilon() );
    }
  }

  
  void create_coordinates_kagome(){
    using namespace std;
    using namespace ffd::user_space;
    using ffd::get, ffd::RealSpace;

    KagomeLattice K({10, 10});

    auto X = CreateCoordinates(K);
    static_assert( is_same<decltype(X),
		   ffd::periodic_coordinate::PeriodicCoordinates<2, int, int, int>>::value );


    check_equality_vectors(RealSpace(X), {0., 0.});

    get<2>(X).Variable = 1;
    check_equality_vectors(RealSpace(X), {1., 0.});

    get<2>(X).Variable = 2;
    check_equality_vectors(RealSpace(X), {.5, .5*sqrt(3.)});

    get<2>(X).Variable = 3; //i.e. = 0
    get<0>(X).Variable = 1;
    check_equality_vectors(RealSpace(X), {2., 0.});

    get<0>(X).Variable = 0;
    get<1>(X).Variable = 1;
    check_equality_vectors(RealSpace(X), {1., sqrt(3.)});

    get<0>(X).Variable = 0;
    get<1>(X).Variable = -1;
    check_equality_vectors(RealSpace(X), {-1., -sqrt(3.)});

    get<0>(X).Variable = 1;
    get<1>(X).Variable = 1;
    get<2>(X).Variable = 1;
    check_equality_vectors(RealSpace(X), {1.+2+1, sqrt(3.)+0+0});

    
    for(int j=-4; j<=5; ++j){
      for(int k=-4; k<=5; ++k){
	for(int u=0; u<3; ++u){
	  get<0>(X).Variable = j;
	  get<1>(X).Variable = k;
	  get<2>(X).Variable = u;
	  check_equality_vectors(RealSpace(X), {2.*j+k+(u==1) + (u==2)*.5, sqrt(3.)*k+(u==2)*.5*sqrt(3.)});
	}
      }
    }

  }

}//namespace
