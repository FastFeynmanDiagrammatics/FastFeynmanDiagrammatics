

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
  NumberParticlesInBlock(int block_number) const{
    auto block_iterator = BlockOracle.Blocks.begin();
    for(int j=0; j < block_number; ++j){
      ++block_iterator;
    }
    
    return std::size(*block_iterator);
  }


}//namespace
