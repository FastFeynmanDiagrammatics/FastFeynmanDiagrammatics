namespace ffd::user_space::imaginary_time_lattice_propagator::unit_test{

  void Nambu_propagator(){

    assert((
	    assertable_Nambu_propagator<3>(
					   {{1, 1, 1}},
					   {{-.5, -.5}},
					   {{
					     {{0., 0., 0.}},
					     {{0., 0., 0.}}
					     }},
					   1.29123,
					   1.4238
					   )
	    ));


    
    assert((
	    assertable_Nambu_propagator<2>(
					   {{2, 3}},
					   {{-.5, -.8}},
					   {{
					     {{1.2121, .81231}},
					     {{.55234, .31231}}
					     }},
					   .89123,
					   3.4238
					   )
	    ));




    
    // std::cerr<<"7bis";
    
    assert((
	    assertable_Nambu_propagator<3>(
					   {{2, 2, 2}},
					   {{-.5, -.8}},
					   {{
					     {{1.2121, .91231, 1.1}},
					     {{.55234, .51231, .8}}
					     }},
					   1.29123,
					   1.4238
					   )
	    ));

  }

}//namespace
