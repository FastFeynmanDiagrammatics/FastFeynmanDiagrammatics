namespace ffd::chebyshev_polynomial_s::unit_test{

  void UnitTest(){

    assert((
	    std::abs(function_approx<30>([](Real x){return 1./(1.+x*x+x*x*x*x);}))
	    < 1e-10
	    ));
    assert((
	    std::abs(function_approx<30>([](Real x){return 1./(1.+x+x*x*x*x);}))
	    < 1e-10
	    ));


    assert(( convolution<40>(2.457)<1e-10 ));

    assert(( convolution<40>(2.157)<1e-10 ));

    assert(( convolution<40, Complex>(1.157)<1e-10 ));
    
  }

}//namespace
