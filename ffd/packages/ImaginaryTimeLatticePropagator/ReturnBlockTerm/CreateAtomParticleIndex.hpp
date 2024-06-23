

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
  CreateAtomParticleIndex(int atom_index, int particle_index) const{
    

    return atom_index +
      number_atoms_unit_cell * particle_index;

    
  }

  
}//namespace
