

namespace ffd::eigenvalues_vectors::unit_test{

  void eigenvalues_vectors_matrix(){
    using ffd::user_space::I;

    
    
    assert((
	    eigenvalues_vectors_matrix_assertable<Real>({{
							  {{1., 2.}},
							  {{2., 5.}}
		}})
	    ));



    
    assert((
	    eigenvalues_vectors_matrix_assertable<Real>({{
							  {{1., -2.}},
							  {{-2., 5.}}
		}})
	    ));


    

    
    assert((
	    eigenvalues_vectors_matrix_assertable<Real>({{
							  {{1., -2.}},
							  {{-2., 1.}}
		}})
	    ));



    assert((
	    eigenvalues_vectors_matrix_assertable<Real>({{
							  {{1., 0.}},
							  {{0., 1.}}
		}})
	    ));



    
    assert((
	    eigenvalues_vectors_matrix_assertable<Real>({{
							  {{-1., 0.}},
							  {{0., 1.}}
		}})
	    ));

    


    
    assert((
	    eigenvalues_vectors_matrix_assertable<Complex>({{
							  {{2., -1.-I}},
							  {{-1.+I, 5.}}
		}})
	    ));



    
    assert((
	    eigenvalues_vectors_matrix_assertable<Complex>({{
							     {{1., 1.-I}},
							     {{1.+I, 1.}}
		}})
	    ));




    assert((
	    eigenvalues_vectors_matrix_assertable<Complex>({{
							     {{0., 1.-I}},
							     {{1.+I, 0.}}
		}})
	    ));



    
    assert((
	    eigenvalues_vectors_matrix_assertable<Complex>({{
							     {{0., -I}},
							     {{I, 0.}}
		}})
	    ));




    assert((
	    eigenvalues_vectors_matrix_assertable<Complex>({{
							     {{1., I}},
							     {{-I, 0.}}
		}})
	    ));



    
    assert((
	    eigenvalues_vectors_matrix_assertable<Complex>({{
							     {{-1., 3.+2.*I}},
							     {{3.-2.*I, -2.}}
		}})
	    ));



    


    
    assert((
	    eigenvalues_vectors_matrix_assertable<Complex>({{
							     {{-1., 3., 2.}},
							     {{3., 2., 4.}},
							     {{2., 4., -3.}}
		}})
	    ));




    assert((
	    eigenvalues_vectors_matrix_assertable<Complex>({{
							     {{1., -1., 5.}},
							     {{-1., 0., -2.}},
							     {{5., -2., -2.}}
		}})
	    ));

    


    Real precision = 1e-8, precision_check = 1e-8;
    assert((
	    eigenvalues_vectors_matrix_assertable<Complex>({{
							     {{0.2      ,-.3-.5*I, -.2+.2*I}},
							     
							     {{-.3+.5*I, .7      , -.8+.4*I}},
							     
							     {{-.2-.2*I   ,-.8-.4*I, -.9}}
		}},
	      precision=12-12,
	      true,
	      precision_check=1e-8)
	    ));



    
    assert((
	    eigenvalues_vectors_matrix_assertable<Real>({{
							     {{1., -1., 5., 4.}},
							     {{-1., 0., -2., -1}},
							     {{5., -2., -2., 2.}},
							     {{4., -1., 2., -1.}}
		}},
	      precision = 1e-12,
	      true,
	      precision_check=1e-8)
	    ));




    precision_check = 1e-5;
    assert((
    	    eigenvalues_vectors_matrix_assertable<Complex>({{
    							     {{.1      , -1.+.2*I, .5+.2*I, .4-.7*I}},
							     
    							     {{-1.-.2*I, -1.     ,-.2+.5*I, -.1-.4*I}},
							     
    							     {{.5-.2*I , -.2-.5*I, -.2    , .3-.2*I }},
							     
    							     {{.4+.7*I , -.1+.4*I, .3+.2*I, -1.     }}
    		}},
    	      precision = 1e-10,
    	      true,
	      precision_check = 1e-8)
    	    ));
    

    
    
    
  }


}//namespace
