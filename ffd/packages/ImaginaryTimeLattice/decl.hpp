namespace ffd::imaginary_time_lattice{

  template<int d, int r, typename... args>
  struct ImaginaryTimeLatticeCoordinatesRecursion{
    using type = typename ImaginaryTimeLatticeCoordinatesRecursion<d, r-1, int, args...>::type;
  };

  template<int d, typename... args>
  struct ImaginaryTimeLatticeCoordinatesRecursion<d, 0, args...>{
    using type = ffd::periodic_coordinate::PeriodicCoordinates<d, Real, args...>;
  };

  


  template<int d, bool has_unit_cell>
  struct ImaginaryTimeLatticeCoordinates;
  

  template<int d>
  struct ImaginaryTimeLatticeCoordinates<d, false>{
    using type = typename ImaginaryTimeLatticeCoordinatesRecursion<d, d>::type;
  };

  
  template<int d>
  struct ImaginaryTimeLatticeCoordinates<d, true>{
    using type = typename ImaginaryTimeLatticeCoordinatesRecursion<d, d+1>::type;
  };
  


}//namespace ffd::imaginary_time_lattice
