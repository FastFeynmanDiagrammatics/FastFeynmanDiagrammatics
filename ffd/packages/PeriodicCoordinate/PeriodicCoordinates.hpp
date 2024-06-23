namespace ffd::periodic_coordinate{

  template<int SpaceDimensions, typename... args>
  using PeriodicCoordinates =
    ffd::class_tuple::ClassTupleTypeInt<PeriodicCoordinate, SpaceDimensions, args...>;
  
}//namespace ffd::periodic_coordinate


