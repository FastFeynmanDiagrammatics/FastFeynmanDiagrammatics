

namespace ffd{

  template<typename T, int SpaceDimensions>
  int sign(ffd::periodic_coordinate::PeriodicCoordinate<T, SpaceDimensions> const& P_){
    if( !P_.LowerUpperBound.has_value() || P_.NotAntiPeriodic ){
      return 1;
    }
    auto const& [Lower, Upper] = P_.LowerUpperBound.value();
    auto const length = Upper - Lower;
    int ret = 1;
    auto variable = P_.Variable;
    while(variable < Lower){
      variable += length;
      ret *= -1;
    }
    while(variable >= Upper){
      variable -= length;
      ret *= -1;
    }
    return ret;
  }


  template<int first, int... last, int SpaceDimensions, typename... args>
  int signSequence(ffd::periodic_coordinate::PeriodicCoordinates<SpaceDimensions, args...> const& P_,
		   ffd::packages_math::StaticRangeSequence<first, first, first, last...> const&){
    return ((ffd::sign(ffd::get<first>(P_)))* ... *(ffd::sign(ffd::get<last>(P_))));
  }

  
  template<int SpaceDimensions, typename... args>
  int sign(ffd::periodic_coordinate::PeriodicCoordinates<SpaceDimensions, args...> const& P_){
    ffd::packages_math::StaticRangeSequence<0, sizeof...(args)> Seq;
    return signSequence(P_, Seq);
  }

  
}//namespace 
