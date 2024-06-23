

namespace ffd::user_space::imaginary_time_lattice_propagator{

  template<int d>
  std::array<Real, d>
  ReturnsWaveVectorFromWaveNumbers(std::array<int, d> const& wave_numbers,
				   std::array<int, d> const& L,
				   std::optional<std::array<Real, d>> const&
				   TwistedBoundaryConditionsVector = {}){
    using ffd::core_math::Pi;
    std::array<Real, d> wave_vector;
    for(int j=0; j < d; ++j){
      wave_vector[j] = 2.*wave_numbers[j]*Pi/L[j];
    }
    if( TwistedBoundaryConditionsVector.has_value() ){
      auto const& twisted_vector = TwistedBoundaryConditionsVector.value();
      for(int j=0; j < d; ++j){
	wave_vector[j] += twisted_vector[j]/L[j];
      }
    }
    return wave_vector;
  }

}//namespace
