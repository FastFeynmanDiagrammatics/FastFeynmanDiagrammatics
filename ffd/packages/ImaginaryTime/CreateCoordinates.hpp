

namespace ffd::user_space{

  template<int space_dimensions=0, bool NotAntiperiodic = true>
  auto CreateCoordinates(ffd::imaginary_time::ImaginaryTime<NotAntiperiodic> const& ImagTime_){
    using ffd::get;
    ffd::periodic_coordinate::PeriodicCoordinates<space_dimensions, Real> tau;
    get<0>(tau).LowerUpperBound = ImagTime_.LowerUpperBound;
    get<0>(tau).Variable = .5*(ImagTime_.LowerUpperBound[0] + ImagTime_.LowerUpperBound[1]);
    get<0>(tau).NotAntiPeriodic = NotAntiperiodic;
    return tau;
  }

}//namespace 
