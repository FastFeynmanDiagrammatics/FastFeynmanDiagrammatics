

namespace ffd::user_space{

  template<typename ImaginaryTimeType,
	   typename LatticeType,
	   typename ActionField,
	   typename PropagatorField>

  Real
  
  ImaginaryTimeLatticePropagator<ImaginaryTimeType,
  				 LatticeType,
  				 ActionField,
  				 PropagatorField>::
  GetInverseTemperatureOfModel() const{
    Real beta = 0;

    
    for( int j: {0, 1} ){
      beta += (2*j-1)*ImagTimeOfModel.LowerUpperBound[j];
    }
    
    
    return beta;
  }

}//namespace
