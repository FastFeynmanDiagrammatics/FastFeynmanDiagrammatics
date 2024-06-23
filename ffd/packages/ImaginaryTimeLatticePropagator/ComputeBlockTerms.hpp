

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
  ComputeBlockTerms(){
    using ffd::vector_range::Range;
    
    for(  int  block:  Range( std::size(BlockOracle.Blocks) )  ){
      auto [block_terms, is_fermion, block_size, not_anomalous] =
	ReturnBlockTerm(block);
      BlocksTerms.push_back(block_terms);
      IsAFermionicBlock.push_back(is_fermion);
      BlockSize.push_back(block_size);
      NotAnomalousBlock.push_back(not_anomalous);
    }
  }

}//namespace
