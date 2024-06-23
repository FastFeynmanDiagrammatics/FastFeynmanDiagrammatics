

namespace ffd::user_space{
    
  template<typename ImaginaryTimeType,
	   typename LatticeType,
	   typename ActionField,
	   typename PropagatorField>
    
  int

  ImaginaryTimeLatticePropagator<ImaginaryTimeType,
				 LatticeType,
				 ActionField,
				 PropagatorField>::
  GetAtomIndex(std::any SpacetimePosition) const{

    
    auto X = std::any_cast<CoordinatesType>(SpacetimePosition);

    
    if constexpr( number_atoms_unit_cell > 1 ){
	return ffd::get<d+1>(X)();
      }else{
      return 0;
    }
    
  }


}//namespace
