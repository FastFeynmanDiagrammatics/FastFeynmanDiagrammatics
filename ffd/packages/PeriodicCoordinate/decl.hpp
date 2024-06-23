namespace ffd::periodic_coordinate{
  
  template<typename VariableType, int space_dimensions>
  class PeriodicCoordinate{
  public:
    VariableType Variable;

    std::optional<std::array<VariableType, 2>> LowerUpperBound;
    
    bool NotAntiPeriodic = true;

    std::vector<std::array<Real, space_dimensions>> RealSpaceVectors;

    std::optional<char> Name;

    bool TranslationInvariant = true;

    PeriodicCoordinate() = default;
    PeriodicCoordinate(const PeriodicCoordinate&) = default;
    PeriodicCoordinate& operator=(const PeriodicCoordinate&) = default;
    
    VariableType operator()() const;
    
  };
  
}//namespace ffd::periodic_coordinate
