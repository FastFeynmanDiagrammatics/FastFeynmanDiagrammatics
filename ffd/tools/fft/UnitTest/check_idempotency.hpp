namespace ffd::fft::unit_test{

  void check_idempotency(){

    assert((
	    idempotency<0>({Complex(1.)})
	    ));

    
    assert((
	    idempotency<1>({Complex(1.), Complex(-1.)})
	    ));



    assert((
	    idempotency<1>({Complex(3.123,1.12321), Complex(1.12312, -1.454)})
	    ));


    assert((
	    idempotency<2>({Complex(2.123,-.12321), Complex(-.2312, 1.54),
				    Complex(-1.34123, -2.1231), Complex(.81276, .912123)})
	    ));


    assert((
	    idempotency<3>({Complex(1.123,-.321), Complex(-.5312, .54),
				    Complex(-1.34123, 2.1231), Complex(.1276, .62123),
				    Complex(1.71, .67), Complex(-2.11, .3),
				    Complex(1.22, .4), Complex(0., 0.)})
	    ));


    assert((
	    idempotency<4>({Complex(1.123,-.321), Complex(-.5312, .54),
				    Complex(-1.34123, 2.1231), Complex(.1276, .62123),
				    Complex(1.71, .67), Complex(-2.11, .3),
				    Complex(1.22, .4), Complex(0., 0.),
				    Complex(-3., 0.), Complex(0., 1.),
				    Complex(.3, .3), Complex(-1.,0.),
				    Complex(0., 0.), Complex(0., -.7)})
	    ));

    

  }

}//namespace
