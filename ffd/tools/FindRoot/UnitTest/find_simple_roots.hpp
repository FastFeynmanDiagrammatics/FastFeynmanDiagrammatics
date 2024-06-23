

namespace ffd::find_root::unit_test{

  void find_simple_roots(){

    Real precision = 1e-12;
    
    auto f = [](Real x) {return x + 2;};

    auto x = FindRealRoot(f,
			  0,
			  precision);

    assert(( x.has_value() ));
    assert((
	    std::abs( x.value() + 2.) < precision
	    ));
    // std::cerr<<x.value()<<std::endl;


    auto f2 = [](Real x) { return std::pow(x+1, 3); };


    std::array<std::optional<Real>, 2> domain;
    domain[0] = 0.;
    assert(( domain[0].has_value() ));
    assert(( !domain[1].has_value() ));
    auto x2 = FindRealRoot(f2,
			   1.,
			   precision,
			   domain);

    
    assert(( !x2.has_value() ));


    
    auto x3 = FindRealRoot(f2,
			   1.,
			   precision);


    assert(( x3.has_value() ));

    
    assert((
	    std::abs(x3.value() + 1.)
	    < precision
	    ));
	    
    
  }

}//namespace
