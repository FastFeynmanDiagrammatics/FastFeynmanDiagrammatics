

namespace ffd::user_space{
  
  template<typename ImaginaryTimeType,
	   typename LatticeType,
	   typename ActionField,
	   typename PropagatorField>
  void
  ImaginaryTimeLatticePropagator<ImaginaryTimeType,
				 LatticeType,		
				 ActionField,
				 PropagatorField>::
  ComputeInverseBravais(){
    std::vector<Real> bravais_vectors(d*d, 0.);
    auto x = CreateCoordinates(LatticeOfModel);
    imaginary_time_lattice_propagator::get_bravais_vectors<d, 0>(x, bravais_vectors);



    InverseBravaisMatrix =
      ffd::inverse_matrix::InverseOfMatrix<Real>(bravais_vectors);
    
  }
  
  
}//namespace
