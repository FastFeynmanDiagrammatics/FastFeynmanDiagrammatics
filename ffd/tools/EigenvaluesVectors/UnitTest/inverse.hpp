

namespace ffd::eigenvalues_vectors::unit_test{

  void inverse(){
    using ffd::user_space::I;

    

    Real precision = 1e-10;
    
    
    assert((
	    inverse_assertable_2x2<Real>({{
				   {{12. , -21.}},
				   {{-21., -2.}}
				   }},
				   precision)
	    ));



    assert((
	    inverse_assertable_2x2<Complex>({{
					   {{12. , -.2+2.*I}},
					   {{-.2-2.*I, -2.}}
		}},
	      precision=1e-10)
	    ));


    

  }
  


}//namespace
