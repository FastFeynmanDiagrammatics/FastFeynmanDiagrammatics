

namespace ffd::user_space{

  template<typename ImaginaryTimeType,
	   typename LatticeType,
	   typename ActionField,
	   typename PropagatorField>
  
  std::vector<int>
  
  ImaginaryTimeLatticePropagator<ImaginaryTimeType,
  				 LatticeType,
  				 ActionField,
  				 PropagatorField>::
  GetLinearSizes() const{
    std::vector<int> L(LatticeType::dimension, 0);

    
    auto const x_linear_sizes = CreateCoordinates(LatticeOfModel);
    imaginary_time_lattice_propagator::
      get_linear_sizes<d, 0>(x_linear_sizes, L);

    
    return L;
  }


 
}//namespace
