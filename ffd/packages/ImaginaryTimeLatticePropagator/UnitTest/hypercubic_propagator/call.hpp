

namespace ffd::user_space::imaginary_time_lattice_propagator::unit_test{

  void hypercubic_propagator(){

    assert((
	    assertable_hypercubic_propagator<3>({{3, 2, 3}},
						{{-.4, -.7}},
						{{
						  {{1.2, 1.5, 1.1}},
						  {{2.11, 1.2, 1.7}}
						  }},
						2.121
						)
	    ));
    
  }


};
