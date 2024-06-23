

namespace ffd::inverse_matrix::unit_test{

  void test_one_minus_A(){
    using ffd::user_space::I;
    
    Real precision = 1e-12;

    
    

    assert((
	    invert_one_minus_A_assertable<Real>({{
						  .02, -.05,
						  .04, -.07
		}},
	      precision)
	    ));




    

    assert((
	    invert_one_minus_A_assertable<Complex>({{
						  +.09+.05*I,-.08+.03*I,+.07-.02*I,
						  
						  +.04+.02*I,+.01-.01*I, 0.       ,
						  
						  -.06+.05*I,+.04-.08*I,-.02+.04*I
		}},
	      precision)
	    ));




    
    assert((
	    invert_one_minus_A_assertable<Complex>({{
						     +.01+.09*I,-.02-.05*I,+.03+.04*I, .09*I,
						  
						     +.07      ,-.05-.02*I,+.03       ,-.04+.03*I,
						  
						     +.03+.07*I,+.02-.06*I,-.01-.03*I,       0.,

						     0.        , .05      ,-.04+.02*I,+.07-.06*I
		}},
	      precision)
	    ));



  }
  

}//namespace
