namespace ffd::fft::unit_test{

  void compare_dft_fft_2d_array(){

    assert((
	    assertable_dft_fft_2d_array<0, 0>({Complex(2.)})
	    ));


    assert((
	    assertable_dft_fft_2d_array<1, 0>({Complex(2., -.3), Complex(-.7, .2)})
	    ));


    assert((
	    assertable_dft_fft_2d_array<0, 1>({Complex(2., -.3), Complex(-.7, .2)})
	    ));


    assert((
	    assertable_dft_fft_2d_array<0, 2>({Complex(2., -.3), Complex(-.7, .2),
					       Complex(.8, -.2), Complex(0., .2)})
	    ));


    assert((
	    assertable_dft_fft_2d_array<2, 0>({Complex(2., -.3), Complex(-.7, .2),
					       Complex(.8, -.2), Complex(0., .2)})
	    ));

    
    assert((
	    assertable_dft_fft_2d_array<1, 1>({Complex(2., -.3), Complex(-.7, .2),
					       Complex(.8, -.2), Complex(0., .2)})
	    ));


    assert((
	    assertable_dft_fft_2d_array<3, 0>({Complex(1.), Complex(3.), Complex(7.), Complex(-2.),
					 Complex(9.), Complex(18.), Complex(-3.222), Complex(1.1222)})
	    ));


    assert((
	    assertable_dft_fft_2d_array<2, 1>({Complex(1.), Complex(3.), Complex(7.), Complex(-2.),
					 Complex(9.), Complex(18.), Complex(-3.222), Complex(1.1222)})
	    ));


    assert((
	    assertable_dft_fft_2d_array<1, 2>({Complex(1.), Complex(3.), Complex(7.), Complex(-2.),
					 Complex(9.), Complex(18.), Complex(-3.222), Complex(1.1222)})
	    ));


    assert((
	    assertable_dft_fft_2d_array<3, 0>({Complex(1.), Complex(3.), Complex(7.), Complex(-2.),
					 Complex(9.), Complex(18.), Complex(-3.222), Complex(1.1222)})
	    ));

    
    assert((
	    assertable_dft_fft_2d_array<4, 0>({Complex(1.), Complex(3.), Complex(7.), Complex(-2.),
					Complex(9.), Complex(18.), Complex(-3.222), Complex(1.1222),
	    Complex(.822), Complex(0.), Complex(1.), Complex(1.2),
	    Complex(9.8), Complex(1.), Complex(0.), Complex(-1.6)})
	    ));


    assert((
	    assertable_dft_fft_2d_array<3, 1>({Complex(1.), Complex(3.), Complex(7.), Complex(-2.),
					Complex(9.), Complex(18.), Complex(-3.222), Complex(1.1222),
	    Complex(.822), Complex(0.), Complex(1.), Complex(1.2),
	    Complex(9.8), Complex(1.), Complex(0.), Complex(-1.6)})
	    ));



    assert((
	    assertable_dft_fft_2d_array<2, 2>({Complex(1.), Complex(3.), Complex(7.), Complex(-2.),
					Complex(9.), Complex(18.), Complex(-3.222), Complex(1.1222),
	    Complex(.822), Complex(0.), Complex(1.), Complex(1.2),
	    Complex(9.8), Complex(1.), Complex(0.), Complex(-1.6)})
	    ));


        assert((
	    assertable_dft_fft_2d_array<1, 3>({Complex(1.), Complex(3.), Complex(7.), Complex(-2.),
					Complex(9.), Complex(18.), Complex(-3.222), Complex(1.1222),
	    Complex(.822), Complex(0.), Complex(1.), Complex(1.2),
	    Complex(9.8), Complex(1.), Complex(0.), Complex(-1.6)})
	    ));


	assert((
	    assertable_dft_fft_2d_array<0, 4>({Complex(1.), Complex(3.), Complex(7.), Complex(-2.),
					Complex(9.), Complex(18.), Complex(-3.222), Complex(1.1222),
	    Complex(.822), Complex(0.), Complex(1.), Complex(1.2),
	    Complex(9.8), Complex(1.), Complex(0.), Complex(-1.6)})
	    ));


	

    
  }


}//namespace
