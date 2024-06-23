

namespace ffd::user_space::imaginary_time_lattice_propagator::unit_test{

  void compute_inverse_bravais(){
    auto G0 = instantiation();

    G0.ComputeInverseBravais();

    for(int j=0; j < 3; ++j){
      for(int k=0; k < 3; ++k){
	assert((
		std::abs(G0.InverseBravaisMatrix[j+k*3] -(j==k)) < 1e-10
		));
      }
    }


    
    
  }


}//namespace
