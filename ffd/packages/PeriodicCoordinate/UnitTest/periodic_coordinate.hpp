

namespace ffd::periodic_coordinate::unit_test{

  void periodic_coordinate(){
    
    PeriodicCoordinate<Real, 3> tau;

    tau.Variable = 2.;

    assert(tau.NotAntiPeriodic);

    tau.LowerUpperBound = {0., 3.};

    tau.Name = 't';

    assert(std::size(tau.RealSpaceVectors) == 0);

    tau.RealSpaceVectors.push_back({.4, .5, .9});
    
    assert( tau.TranslationInvariant );
    
  }

}//namespace
