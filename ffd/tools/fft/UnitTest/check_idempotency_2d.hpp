namespace ffd::fft::unit_test{

  void check_idempotency_2d(){

    assert((
	    idempotency_2d<0, 0>({Complex(2., .2)})
	    ));

    assert((
	    idempotency_2d<1, 0>({Complex(1., -.2), Complex(-1., .9)})
	    ));


    assert((
	    idempotency_2d<0, 1>({Complex(1.), Complex(-1., .4)})
	    ));


    assert((
	    idempotency_2d<2, 0>({Complex(1., -.2), Complex(-1., .9),
			       Complex(.3, -.4), Complex(-.8, .2)})
	    ));

    assert((
	    idempotency_2d<1, 1>({Complex(1., -.2), Complex(-1., .9),
			       Complex(.3, -.4), Complex(-.8, .2)})
	    ));

    assert((
	    idempotency_2d<0, 2>({Complex(1., -.2), Complex(-1., .9),
			       Complex(.3, -.4), Complex(-.8, .2)})
	    ));


    assert((
	    idempotency_2d<3, 0>({Complex(1., -.2), Complex(-1., .9),
			       Complex(.3, -.4), Complex(-.8, .2),
			       Complex(0.), Complex(0., 2.3),
			       Complex(-.3, .7), Complex(-.8, .1)})
	    ));

        assert((
	    idempotency_2d<2, 1>({Complex(1., -.2), Complex(-1., .9),
			       Complex(.3, -.4), Complex(-.8, .2),
			       Complex(0.), Complex(0., 2.3),
			       Complex(-.3, .7), Complex(-.8, .1)})
	    ));

	assert((
	    idempotency_2d<1, 2>({Complex(1., -.2), Complex(-1., .9),
			       Complex(.3, -.4), Complex(-.8, .2),
			       Complex(0.), Complex(0., 2.3),
			       Complex(-.3, .7), Complex(-.8, .1)})
	    ));
    
    assert((
	    idempotency_2d<0, 3>({Complex(1., -.2), Complex(-1., .9),
			       Complex(.3, -.4), Complex(-.8, .2),
			       Complex(0.), Complex(0., 2.3),
			       Complex(-.3, .7), Complex(-.8, .1)})
	    ));


        assert((
		idempotency_2d<4, 0>({Complex(1., -.2), Complex(-1., .9),
			       Complex(.3, -.4), Complex(-.8, .2),
			       Complex(0.), Complex(0., 2.3),
			       Complex(-.3, .7), Complex(-.8, .1),
			       
			       Complex(.3, .8), Complex(-.2, -.2),
			       Complex(.7, .8), Complex(0.),
			       Complex(.4, .2), Complex(-.1, 0.3),
			       Complex(-.7, -.8), Complex(.3, 0.)})
	    ));
        assert((
		idempotency_2d<3, 1>({Complex(1., -.2), Complex(-1., .9),
			       Complex(.3, -.4), Complex(-.8, .2),
			       Complex(0.), Complex(0., 2.3),
			       Complex(-.3, .7), Complex(-.8, .1),
			       
			       Complex(.3, .8), Complex(-.2, -.2),
			       Complex(.7, .8), Complex(0.),
			       Complex(.4, .2), Complex(-.1, 0.3),
			       Complex(-.7, -.8), Complex(.3, 0.)})
	    ));
        assert((
		idempotency_2d<2, 2>({Complex(1., -.2), Complex(-1., .9),
			       Complex(.3, -.4), Complex(-.8, .2),
			       Complex(0.), Complex(0., 2.3),
			       Complex(-.3, .7), Complex(-.8, .1),
			       
			       Complex(.3, .8), Complex(-.2, -.2),
			       Complex(.7, .8), Complex(0.),
			       Complex(.4, .2), Complex(-.1, 0.3),
			       Complex(-.7, -.8), Complex(.3, 0.)})
	    ));
        assert((
		idempotency_2d<1, 3>({Complex(1., -.2), Complex(-1., .9),
			       Complex(.3, -.4), Complex(-.8, .2),
			       Complex(0.), Complex(0., 2.3),
			       Complex(-.3, .7), Complex(-.8, .1),
			       
			       Complex(.3, .8), Complex(-.2, -.2),
			       Complex(.7, .8), Complex(0.),
			       Complex(.4, .2), Complex(-.1, 0.3),
			       Complex(-.7, -.8), Complex(.3, 0.)})
	    ));
        assert((
		idempotency_2d<0, 4>({Complex(1., -.2), Complex(-1., .9),
			       Complex(.3, -.4), Complex(-.8, .2),
			       Complex(0.), Complex(0., 2.3),
			       Complex(-.3, .7), Complex(-.8, .1),
			       
			       Complex(.3, .8), Complex(-.2, -.2),
			       Complex(.7, .8), Complex(0.),
			       Complex(.4, .2), Complex(-.1, 0.3),
			       Complex(-.7, -.8), Complex(.3, 0.)})
	    ));



  }
  
}//namespace
