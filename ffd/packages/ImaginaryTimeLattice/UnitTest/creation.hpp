

namespace ffd::imaginary_time_lattice::unit_test{

  void creation(){
    using namespace ffd::user_space;
    using ffd::get, ffd::sign;
    
    ffd::imaginary_time::AntiperiodicImaginaryTime ImagTime(2.);

    auto tau = ffd::user_space::CreateCoordinates<2>(ImagTime);

    ffd::lattice::HypercubicLattice<2> square({10, 10});

    auto x = CreateCoordinates(square);

    auto X = tau*x;


    ffd::lattice::KagomeLattice K({10, 10});

    auto y = CreateCoordinates(K);

    auto Y = tau*y;

    assert( std::abs( get<0>(Y)() - 1. ) < std::numeric_limits<Real>::epsilon() );
    assert( get<1>(Y)() == 0 );
    assert( get<2>(Y)() == 0 );
    assert( get<3>(Y)() == 0 );

    auto Sign = sign(Y);
    assert(Sign==1);
    auto r = RealSpace(Y);   
    for(int j=0; j < 2; ++j){
      assert( std::abs(r[j]) < std::numeric_limits<Real>::epsilon());
    }


    
    using coordinates = typename ImaginaryTimeLatticeCoordinates<3, true>::type;

    coordinates Z;

    static_assert( std::is_same<coordinates,
		   ffd::periodic_coordinate::PeriodicCoordinates<3,
		   Real,
		   int,
		   int,
		   int,
		   int>>::value );


    ImaginaryTimeLatticeCoordinates<3, true>::type Z2;

    static_assert( std::is_same<decltype(Z),
		   decltype(Z2)>::value );

    

    
    
    
  }

}//namespace
