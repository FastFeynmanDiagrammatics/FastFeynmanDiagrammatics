

namespace ffd::fft::unit_test{

  void
  check_invert_binary_number(){

    assert((
	    0 == invert_binary_number(0u, 2)
	    ));


    assert((
	    0 == invert_binary_number(0u, 1)
	    ));


    
    assert((
	    0 == invert_binary_number(0u, 4)
	    ));


    assert((
	    1 == invert_binary_number(1, 1)
	    ));

    
    assert((
	    2 == invert_binary_number(1, 2)
	    ));

    
    assert((
	    10 == invert_binary_number(5, 4)
	    ));

    
    assert((
	    3 == invert_binary_number(12, 4)
	    ));


    
    int const n_digits = 5;
    for(int j=0; j<n_digits; ++j){
      assert((
	      1<<(n_digits-j-1) == invert_binary_number(1<<j, n_digits)
	      ));
    }


    
    for(int j=0; j< (1<<n_digits); ++j){
      auto jj = invert_binary_number( invert_binary_number(j, n_digits), n_digits );
      assert((
	      j == jj
	      ));
    }


    
    
  }


}//namespace
