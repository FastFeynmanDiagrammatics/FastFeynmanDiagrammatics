

namespace ffd::periodic_coordinate::unit_test{

  void tau_real_space(){
    using ffd::RealSpace;
    
    PeriodicCoordinate<Real, 4> tau;

    for(int j=0; j<4; ++j){
      assert(std::abs(RealSpace(tau)[j]) < 10*std::numeric_limits<Real>::epsilon() );
    }

  }

}//namespace
