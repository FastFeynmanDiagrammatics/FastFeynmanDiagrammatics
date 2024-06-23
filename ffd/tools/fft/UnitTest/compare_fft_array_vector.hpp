namespace ffd::fft::unit_test{

  void compare_fft_array_vector(){

    assert((
	    assertable_dft_fft_array<0>({Complex(1.)})
	    ));


    assert((
	    assertable_dft_fft_array<1>({Complex(1.), Complex(3.,-1.)})
	    ));
    


    assert((
	    assertable_dft_fft_array<2>({Complex(1.), Complex(3.), Complex(7.), Complex(-2.)})
	    ));


    assert((
	    assertable_dft_fft_array<3>({Complex(1.), Complex(3.), Complex(7.), Complex(-2.),
					 Complex(9.), Complex(18.), Complex(-3.222), Complex(1.1222)})
	    ));


    assert((
	    assertable_dft_fft_array<4>({Complex(1.), Complex(3.), Complex(7.), Complex(-2.),
					Complex(9.), Complex(18.), Complex(-3.222), Complex(1.1222),
	    Complex(.822), Complex(0.), Complex(1.), Complex(1.2),
	    Complex(9.8), Complex(1.), Complex(0.), Complex(-1.6)})
	    ));


  }

}//namespace
