

namespace ffd::fft::unit_test{

  void compare_dft_fft_2d(){

    
    assert((
	    assertable_dft_fft_2d({Complex(1., 0.)}, 1) 
	    ));

    

    assert((
	    assertable_dft_fft_2d({Complex(1., 0.3), Complex(2., 0.5)}, 1) 
	    ));


    
    assert((
	    assertable_dft_fft_2d({Complex(1., 0.3), Complex(2., 0.5)}, 2) 
	    ));


    
    assert((
	    assertable_dft_fft_2d({Complex(1., 0.3), Complex(2., 0.5),
				   Complex(.2, .5), Complex(.1, -.4)}) 
	    ));


    
    assert((
	    assertable_dft_fft_2d({Complex(1., 0.3), Complex(2., 0.5),
				   Complex(.2, .5), Complex(.1, -.4)}, 1) 
	    ));


    
    assert((
	    assertable_dft_fft_2d({Complex(1., 0.3), Complex(2., 0.5),
				   Complex(.2, .5), Complex(.1, -.4)}, 4) 
	    ));


    
    assert((
	    assertable_dft_fft_2d({Complex(1., 0.3), Complex(2., 0.5),
				   Complex(.2, .5), Complex(.1, -.4),
				   Complex(-.5, .9), Complex(-1, .4),
				   Complex(.2, .5), Complex(.1, -.4)}, 2)
	    ));



    
    assert((
    	    assertable_dft_fft_2d({Complex(1., 0.3), Complex(2., 0.5),
    				   Complex(.2, .5), Complex(.1, -.4),
				   
    				   Complex(-.5, .9), Complex(-1, .4),
    				   Complex(.2, .5), Complex(.1, -.4),
				   
    				   Complex(.3, .4), Complex(-.1, 0),
    				   Complex(.2, .8), Complex(-.2, -.4),
				   
    				   Complex(.1, .3), Complex(.1, .1),
    				   Complex(.9, -.3), Complex(.3, .5)}, 2)
    	    ));
    

    
    assert((
    	    assertable_dft_fft_2d({Complex(1., 0.3), Complex(2., 0.5),
    				   Complex(.2, .5), Complex(.1, -.4),
				   
    				   Complex(-.5, .9), Complex(-1, .4),
    				   Complex(.2, .5), Complex(.1, -.4),
				   
    				   Complex(.3, .4), Complex(-.1, 0),
    				   Complex(.2, .8), Complex(-.2, -.4),
				   
    				   Complex(.1, .3), Complex(.1, .1),
    				   Complex(.9, -.3), Complex(.3, .5)})
    	    ));


    
        assert((
	    assertable_dft_fft_2d({Complex(1., 0.3), Complex(2., 0.5),
				   Complex(.2, .5), Complex(.1, -.4),
				   
				   Complex(-.5, .9), Complex(-1, .4),
				   Complex(.2, .5), Complex(.1, -.4),
				   
				   Complex(.3, .4), Complex(-.1, 0),
				   Complex(.2, .8), Complex(-.2, -.4),
				   
				   Complex(.1, .3), Complex(.1, .1),
				   Complex(.9, -.3), Complex(.3, .5)}, 1)
	    ));




	assert((
	    assertable_dft_fft_2d({Complex(1., 0.3), Complex(2., 0.5),
				   Complex(.2, .5), Complex(.1, -.4),
				   
				   Complex(-.5, .9), Complex(-1, .4),
				   Complex(.2, .5), Complex(.1, -.4),
				   
				   Complex(.3, .4), Complex(-.1, 0),
				   Complex(.2, .8), Complex(-.2, -.4),
				   
				   Complex(.1, .3), Complex(.1, .1),
				   Complex(.9, -.3), Complex(.3, .5)}, 16)
	    ));


	
  }


}//namespace
