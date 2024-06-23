

namespace ffd::user_space{

  template<typename ImaginaryTimeType,
	   typename LatticeType,
	   typename ActionField,
	   typename PropagatorField>
    
  auto

  ImaginaryTimeLatticePropagator<ImaginaryTimeType,
				 LatticeType,
				 ActionField,
				 PropagatorField>::
  GetRealSpacePositions(std::any SpacetimePosition) const{
    using ffd::vector_range::Range;

      
    std::array<Real, d> real_positions;
    real_positions.fill(0.);
      
    auto X = std::any_cast<CoordinatesType>(SpacetimePosition);
    auto R = RealSpace(X);
    for(int u: Range(d)){
      for(int v: Range(d)){
	real_positions[u] += InverseBravaisMatrix[u*d+v]*R[v];
      }
    }

    return real_positions;
  }


}//namespace
