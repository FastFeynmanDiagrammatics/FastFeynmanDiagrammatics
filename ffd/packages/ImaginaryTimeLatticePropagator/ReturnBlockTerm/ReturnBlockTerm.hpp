

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
  ReturnBlockTerm(int block_number) const{
    using ffd::vector_range::Range;

    
    std::vector<BlockTerm> 
      value_positions_indexes_isdaggers;

    
    std::optional<bool> is_fermionic_optional;

    
    int block_size =
      NumberParticlesInBlock(block_number)*number_atoms_unit_cell;
    

    std::optional<bool> NotAnomalous;
    

    for(auto const& bilinear: Action){
      if(  BlockOracle( bilinear[0][0] )  ==  block_number  ){

	
	auto scalar_value = bilinear.Value;

	
	std::array<std::array<Real, d>, 2> realspace_positions;
	for(int dot: {0, 1}){
	  realspace_positions[dot] = GetRealSpacePositions( bilinear[dot].Position );
	}
	
	
	std::array<int, 2> particle_indexes;
	for(int dot: {0, 1}){
	  particle_indexes[dot] = GetParticleIndex( bilinear[dot][0] );
	}


	std::array<int, 2> atom_indexes;
	for(int dot: {0, 1}){
	  atom_indexes[dot] = GetAtomIndex( bilinear[dot].Position );
	}


	std::array<int, 2> indexes;
	for(int dot: {0, 1}){
	  indexes[dot] =
	    CreateAtomParticleIndex(atom_indexes[dot],
				    particle_indexes[dot]);
	}

	
	if( !is_fermionic_optional.has_value() ){
	  is_fermionic_optional = bilinear[0][0].IsFermion();
	}
	
	
	std::array<std::vector<bool>, 2> is_dagger_sets;
	for(int dot: {0, 1}){
	  
	  is_dagger_sets[dot] = imaginary_time_lattice_propagator::
	    IsDagger_set( bilinear[dot][0].Dagger() );
	  
	  scalar_value  /=  ( std::size( is_dagger_sets[dot] ) == 2 ? std::sqrt(2) : 1. );
	  
	}

	
	for( auto is_daggers: ffd::cartesian_product::CartesianProduct(is_dagger_sets) ){

	  auto element =
	    std::tuple{scalar_value,
		       realspace_positions,
		       indexes,
		       is_daggers};
	    
	  value_positions_indexes_isdaggers.push_back(element);

	    
	  if(is_daggers[0] == is_daggers[1] && !NotAnomalous.has_value()){
	    NotAnomalous = false;
	  }
	  
	}
      }
    }

    return std::tuple{value_positions_indexes_isdaggers,
	is_fermionic_optional.value(),
	block_size,
	NotAnomalous};
  }

}//namespace
