

namespace ffd::imaginary_time::unit_test{

  void create_coordinates(){
    using ffd::get;
    PeriodicImaginaryTime It(2.);
    auto tau = ffd::user_space::CreateCoordinates(It);

    assert( std::abs( get<0>(tau)() - 1 ) < std::numeric_limits<Real>::epsilon() );

    for(int j=0; j < 2; ++j){
      assert( std::abs( get<0>(tau).LowerUpperBound.value()[j] - (j==1)*2) <
	      std::numeric_limits<Real>::epsilon() );
    }
    
  }

}//namespace
