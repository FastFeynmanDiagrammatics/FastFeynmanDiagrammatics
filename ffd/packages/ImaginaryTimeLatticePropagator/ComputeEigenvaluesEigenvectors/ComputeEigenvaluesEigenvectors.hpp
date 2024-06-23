

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
  ComputeEigenvaluesEigenvectors(){
    using ffd::vector_range::Range;

    
    auto L = GetLinearSizes();
    std::array<std::vector<Real>, d> wave_vector_ranges;
    for( int j: Range(d) ){
      wave_vector_ranges[j] = imaginary_time_lattice_propagator::
	ReturnWaveVectorRange(L[j], TwistedBoundaryConditionsVector[j]);
    }

    
    int const number_of_blocks = std::size(BlockOracle.Blocks);
    EigenvaluesVectorsOfBlock1Momentum2.resize( number_of_blocks );

    
    for( int block_number: Range(number_of_blocks) ){
      auto const block_term =          BlocksTerms[block_number];
      auto const is_fermion =    IsAFermionicBlock[block_number];
      auto const block_size =            BlockSize[block_number];
      auto const not_anomalous = NotAnomalousBlock[block_number];
      auto & eigenvalues_vec = EigenvaluesVectorsOfBlock1Momentum2[block_number];
	
      for( auto wave_vector: ffd::cartesian_product::
	     CartesianProduct(wave_vector_ranges) ){

	auto H_k = ReturnBlockHamiltonian_k(block_term,
					    is_fermion,
					    block_size,
					    wave_vector);
	
	
	if( !not_anomalous.has_value() ){
	  H_k = imaginary_time_lattice_propagator::HalveMatrix(H_k);
	}
	

	eigenvalues_vec.push_back( ffd::eigenvalues_vectors::
				   EigenvaluesVectors<Complex>(H_k) );
	
      }
    }
  }


}//namespace
