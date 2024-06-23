

namespace ffd::user_space{

  template<typename ImaginaryTimeType,
	   typename LatticeType,
	   typename ActionField,
	   typename PropagatorField>
  
  PropagatorField
  
  ImaginaryTimeLatticePropagator<ImaginaryTimeType,
  				 LatticeType,
  				 ActionField,
  				 PropagatorField>::
  operator()(ffd::feynman_edge::
	     QuantumFieldPositions Psi_X) const{
    using ffd::quantum_field::QuantumField;
    using ffd::vector_range::Range;


    bool const is_fermion = std::get<QuantumField>( Psi_X[0] ).IsFermion();
    Real const zeta  =  ( is_fermion ? -1 : 1 );

    
    std::array<std::vector<bool>, 2> is_daggers_sets;
    for( int dot: {0, 1} ){
      is_daggers_sets[dot] = imaginary_time_lattice_propagator::
	IsDagger_set( std::get<QuantumField>( Psi_X[dot] ).Dagger() );
    }


    if(  is_daggers_sets[0]  >  is_daggers_sets[1]  ){
      return zeta*(*this)(  imaginary_time_lattice_propagator::
			    SwapArrayElements( Psi_X )  );
    }

    std::array<int, 2> particle_indexes;
    for( int dot: {0, 1} ){
      particle_indexes[dot] =
	GetParticleIndex(  std::get<QuantumField>( Psi_X[dot] )  );
    }


    std::array<int, 2> atom_indexes;
    for( int dot: {0, 1} ){
      atom_indexes[dot] =
	GetAtomIndex(  std::get<std::any>( Psi_X[dot] )  );
    }
    
    
    std::array<int, 2> indexes;
    for( int dot: {0, 1} ){
      indexes[dot] =
	CreateAtomParticleIndex(atom_indexes[dot],
				particle_indexes[dot]);
    }
    
    
    std::array<int, 2> block_numbers;
    for( int dot: {0, 1} ){
      block_numbers[dot] = BlockOracle(  std::get<QuantumField>( Psi_X[dot] )  );
    }
    int const block_number = (block_numbers[0] == block_numbers[1] ?
			      block_numbers[0] : -1);
    if(block_number == -1){
      return 0.;
    }

    
    auto delta_R = imaginary_time_lattice_propagator::
      RealSpaceDistance<d, CoordinatesType>( Psi_X );
    decltype(delta_R) delta_r;
    delta_r.fill(0.);
    for( int j: Range(d) ){
      for( int m: Range(d) ){
	delta_r[j] += InverseBravaisMatrix[j*d + m]*delta_R[m];
      }
    }

    
    auto delta_tau = imaginary_time_lattice_propagator::
      ImaginaryTimeDistance<CoordinatesType>( Psi_X );


    Real Beta = GetInverseTemperatureOfModel();

    
    Complex value_propagator = imaginary_time_lattice_propagator::
      Constrain1BetweenZeroAnd2ReturnSign(delta_tau, Beta, is_fermion);

    



    
    
    if( !NotAnomalousBlock[block_number].has_value() ){
      if( is_daggers_sets[0] ==
	  is_daggers_sets[1] ){
	return 0.;
      }
    }


    auto const& eigen_val_vec_block =
      EigenvaluesVectorsOfBlock1Momentum2[block_number];
    
    

    auto L = GetLinearSizes();
    std::array<std::vector<Real>, d> wave_vector_ranges;
    for( int j: Range(d) ){
      wave_vector_ranges[j] = imaginary_time_lattice_propagator::
	ReturnWaveVectorRange(L[j], TwistedBoundaryConditionsVector[j]);
    }

    Complex sum_over_k = 0;
    int iterator_k = 0;
    for( auto wave_vector: ffd::cartesian_product::
	   CartesianProduct(wave_vector_ranges) ){
      
      
      auto const& [eigenvalues_k,
		   eigenvectors_k] =
	eigen_val_vec_block[iterator_k];
      ++iterator_k;
	

      Complex expI_phase_k =
	expI( imaginary_time_lattice_propagator::
	      ScalarProductArrays<d>(wave_vector, delta_r) );

      for( auto u: Range( std::size(eigenvalues_k) ) ){

	
	Real g_value = imaginary_time_lattice_propagator::
	  GreenFunctionZetaEnergyTau(is_fermion,
				     Beta,
				     eigenvalues_k[u],
				     delta_tau);
	
	
	for( auto is_daggers: ffd::cartesian_product::
	       CartesianProduct(is_daggers_sets) ){

	
	  std::array<int, 2> matrix_indexes;
	  for( int dot: {0, 1} ){
	    matrix_indexes[dot] = indexes[dot];
	    if( is_daggers[dot] == 1-dot ){
	      matrix_indexes[dot] += BlockSize[block_number];
	      assert(( std::size_t(matrix_indexes[dot]) < size(eigenvectors_k) ));
	    }
	  }


	  Complex product_vectors = 1.;
	  for( int dot: {0, 1}){
	    Complex const& component_eigenvector = eigenvectors_k[u][ matrix_indexes[dot] ];
	    if(dot == 0){
	      product_vectors *= component_eigenvector;
	    }else{
	      product_vectors *= std::conj( component_eigenvector );
	    }
	  }


	  for( int dot: {0, 1} ){
	    if( std::size(is_daggers_sets[dot]) == 2){
	      product_vectors /= std::sqrt(2);
	    }
	  }

	  
	  sum_over_k += expI_phase_k*g_value*product_vectors;

	  
	}
      }
    }
    

    value_propagator *= .5*sum_over_k;
    for( int j: Range(d) ){
      value_propagator /= L[j];
    }


    if constexpr( std::is_same_v<PropagatorField, Real> ){
return std::real(value_propagator);
}else{
      return value_propagator;
    }
  }

}//namespace
