

namespace ffd::inverse_matrix::unit_test{

  void inverse_2x2(){
    using ffd::user_space::I;
    
    Real precision = 1e-12;
    

    assert((
	    inverse_2x2_assertable<Real>({{
				     1., 2.,
				     3., 4.
		}})
	    ));

    


    assert((
	    inverse_2x2_assertable<Real>({{
				     .12311312, .0,
				     .5823213, .002
		}},
	      precision)
	    ));



    

    assert((
	    inverse_2x2_assertable<Complex>({{
				     -.4123978+.12378*I, .324789-.768132*I,
				     -.823213-.121*I, .128398-.0123*I
		}},
	      precision)
	    ));


    
    
				     


  }


}//Namespace
