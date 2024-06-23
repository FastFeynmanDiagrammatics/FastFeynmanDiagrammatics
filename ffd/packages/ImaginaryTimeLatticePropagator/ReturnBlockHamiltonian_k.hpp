

namespace ffd::user_space{

  template<typename ImaginaryTimeType,
	   typename LatticeType,
	   typename ActionField,
	   typename PropagatorField>
  
  std::vector<std::vector<Complex>>
  
  ImaginaryTimeLatticePropagator<ImaginaryTimeType,
  				 LatticeType,
  				 ActionField,
  				 PropagatorField>::
  ReturnBlockHamiltonian_k(std::vector<
			   std::tuple<			      
			   ActionField,
			   std::array<std::array<Real, d>, 2>,
			   std::array<int, 2>,
			   std::array<bool, 2>>>
			   const&
			   terms_quadratic_action_block,
			   bool is_fermionic,
			   int block_size,
			   std::array<Real, d>
			   wave_vector) const{
    using ffd::vector_range::Range;
    std::size_t const linear_matrix_size = block_size*2;
    
    
    std::vector<std::vector<Complex>>
      hamiltonian(linear_matrix_size,
		  std::vector<Complex>(linear_matrix_size, 0.));
    
    
    for( auto const& bilinear: terms_quadratic_action_block ){
      auto const& [t, r, indexes, b] = bilinear;
      
      
      std::array<Real, d> r0_r1;
      r0_r1.fill(0.);
      for(int j: {0, 1}){
	for(int m: Range(d)){ 
	  r0_r1[m] += (1-2*j)*r[j][m];
	}
      }
      

      Real k__dot__r0_r1 = 0;
      for( int m: Range(d) ){
	k__dot__r0_r1 += r0_r1[m]*wave_vector[m];
      }
      
      
      for(int j: {0, 1}){
	int sign  =  (j == 0  ?  1  :  1 - 2*is_fermionic);
	
	hamiltonian
	  [  indexes[j]    +  block_size * ( 1 - b[j] )  ]
	  [  indexes[1-j]  +  block_size * b[1-j]      ]
	  +=  t*(1.*sign)*ffd::user_space::expI(-(1-2*j)*k__dot__r0_r1);
	
      }
    }
    return hamiltonian; 
  }
  
}//namespace
