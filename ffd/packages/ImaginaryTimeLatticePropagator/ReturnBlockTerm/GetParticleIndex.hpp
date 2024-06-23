

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
  GetParticleIndex(ffd::quantum_field::QuantumField field) const{
    int block_number = BlockOracle(field);


    auto block_iterator = BlockOracle.Blocks.begin();
    for(int j=0; j < block_number; ++j){
      ++block_iterator;
    }

    
    auto unbar_field  =  field.Dagger() >= 0  ?
      field : Bar(field);

    
    return ffd::relation_sets::
      PositionOfElement1InSet2(unbar_field, *block_iterator).value();
  }


}//namespace
