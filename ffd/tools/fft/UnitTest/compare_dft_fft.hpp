namespace ffd::fft::unit_test{

  void
  compare_dft_fft(){
    
    assert((
	    assertable_dft_fft({1.})
	    ));


    assert((
	    assertable_dft_fft({3.})
	    ));


    
    assert((
	    assertable_dft_fft({{Complex(3.), Complex(4.)}})
	    ));

    

    assert((
	    assertable_dft_fft({{Complex(1., .2), Complex(-1., .44), Complex(4., 2.), Complex(.6, .5)}})
	    ));


    assert((
	    assertable_dft_fft({{Complex(1., .2), Complex(-1., .44), Complex(4., 2.), Complex(.6, .5),
				 Complex(0.7, .9), Complex(-0.2, .2), Complex(-1., 1.), Complex(2.2, .12)
		}})
	    ));

    

    
    {
      int const log2_size = 4;
      unsigned long const size_dft = 1<<log2_size;
      std::vector<Complex> vec(size_dft);
      using ffd::random_distributions::Proba;
      for( int j=0; j< size_dft; ++j){
	vec[j] = Complex(2*Proba()-1, 2*Proba()-1);
      }
      assert((
	      assertable_dft_fft(vec)
	      ));
    }

    

    
    {
      int const log2_size = 6;
      unsigned long const size_dft = 1<<log2_size;
      std::vector<Complex> vec(size_dft);
      using ffd::random_distributions::Proba;
      for( int j=0; j< size_dft; ++j){
	vec[j] = Complex(2*Proba()-1, 2*Proba()-1);
      }
      assert((
	    assertable_dft_fft(vec)
	    ));
    }
    

    
  }
    
}//namespace
