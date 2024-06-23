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
  BuildChebyshev(){
    using ffd::vector_range::Range;
    using ffd::quantum_field::QuantumField;

    
    for(  int block_number:
	    Range( size(BlockOracle.Blocks) )  ){
      auto iter_blocks = BlockOracle.Blocks.begin();
      std::advance(iter_blocks, block_number);
      auto const& block = *iter_blocks;
      
      
      for(   auto nondagger_field_pair:
	       CartesianProduct<2>(  block  )   ){
	
	
	std::array<std::set<QuantumField>, 2> field_sets;
	std::array<bool, 2> second_part_of_matrix;
	for( int j: {0, 1} ){
	  field_sets[j].insert(
			       j==0  ?
			       nondagger_field_pair[0]  :
			       Bar( nondagger_field_pair[1] )
			       );
	  second_part_of_matrix[j] = false;
	  if( NotAnomalousBlock[block_number].has_value() &&
	      nondagger_field_pair[j].Dagger() != 0 ){
	    second_part_of_matrix[j] = true;
	    field_sets[j].insert(
				 j==0  ?
				 Bar( nondagger_field_pair[0] )  :
				 nondagger_field_pair[1]
				 );
	  }
	}
	
	
	for( auto left_right_fields:
	       CartesianProduct(field_sets) ){

	  
	  std::array<int, 2> particle_indexes;
	  for( int j: Range(2) ){
	    particle_indexes[j]  =
	      GetParticleIndex( left_right_fields[j] );
	  }
	  
	  
	  std::array<int, d> L;
	  auto L_vec = GetLinearSizes();
	  for( int j: Range(d) ){  L[j] = L_vec[j];  }
	  

	  auto O =  CreateCoordinates(ImagTimeOfModel, LatticeOfModel);
	  component<0>(O) = std::numeric_limits<Real>::min();
	  
	  
	  std::array<decltype(O), 2> X;
	  for( int left_right: Range(2) ){
	    X[left_right] = O;
	  }
	  
	  
	  for(  auto left_right_atoms:
		  CartesianProductRange<2>( number_atoms_unit_cell )  ){
	    

	    std::array<int, 2> indexes;
	    for( int j: Range(2) ){
	      indexes[j] =
		CreateAtomParticleIndex(left_right_atoms[j],
					particle_indexes[j]);
	    }
	    
	    
	    if constexpr(number_atoms_unit_cell > 1){
		for( int left_right: Range(2) ){
		  component<d+1>(X[left_right]) = left_right_atoms[left_right];
		}
	      }

	    
	    for(  auto r: CartesianProductRange( L )  ){
	      

	      imaginary_time_lattice_propagator::
		SetSpaceCoordinates<d>(X[0], r);
	      
	      
	      std::array<std::pair<QuantumField, std::any>, 2> left_right_field_X;
	      for( int left_right: Range(2) ){
		left_right_field_X[left_right] =
		  std::pair{ left_right_fields[left_right],
			     std::any(X[left_right]) };
	      }
	      auto X_tau = X[0];
	      
	      
	      auto G0_tau = [this,
			     XX = X_tau,
			     &left_right_field_X](Real tau) mutable ->PropagatorField{
			      component<0>(XX) += tau;
			      std::get<1>(left_right_field_X[0]) = XX;
			      return (*this)(left_right_field_X);
			    };
	      
	      
	      Real Beta = 
		GetInverseTemperatureOfModel();
	      auto cheby_poly = ffd::chebyshev_polynomial::
		ChebyshevPolynomial<PropagatorField>(G0_tau,
						     {0, Beta},
						     AbsolutePrecision);
	      

	      ChebyshevPropagator.insert_or_assign(std::tuple{
		  left_right_fields,
		    left_right_atoms,
		    r},
		cheby_poly);
		    
	      
	    }

	  }
	  
	}
	       
      }
    }
  }





}//namespace
