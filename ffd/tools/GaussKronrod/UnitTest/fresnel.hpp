namespace ffd::gauss_kronrod::unit_test{

  void integrate_fresnel(){

    Real precision = 1e-12;

    std::array<Real, 1ul<<10> big_array;
    big_array.fill(2.);

    std::array<Real, 2> interval{-1., 4.};

    
    auto fresnel = [=](Real x){ return std::sin(big_array[1ul<<5]*x*x);};

    
    // long double const exact_fresnel =
    //   0.535613733664510686965298555078591800814319130364619098357l;
    
    long double const exact_fresnel = 0.889287667403153461014670771020019798374872287683738414360l;

    assert((
    	    std::abs(
    		     Integrate(fresnel, interval, precision)-
    		     exact_fresnel
    		     ) < precision
	    
    	    ));


    // std::cerr<<"exact    = "<<std::setprecision(16)<<exact_fresnel<<std::endl;
        
    
    // std::cerr<<"adaptive = "<<std::setprecision(16)<<Integrate(fresnel, interval)<<std::endl;

    // auto [val3, err3] = GaussKronrodIntegral<3>(fresnel, interval);
    // std::cerr<<"GK3      = "<<std::setprecision(16)<<val3<<" "<<err3<<std::endl;
    

    // auto [val7, err7] = GaussKronrodIntegral<7>(fresnel, interval);
    // std::cerr<<"GK7      = "<<std::setprecision(16)<<val7<<" "<<err7<<std::endl;


    // auto [val15, err15] = GaussKronrodIntegral<15>(fresnel, interval);
    // std::cerr<<"GK15     = "<<std::setprecision(16)<<val15<<" "<<err15<<std::endl;


    // auto [val31, err31] = GaussKronrodIntegral<31>(fresnel, interval);
    // std::cerr<<"GK31     = "<<std::setprecision(16)<<val31<<" "<<err31<<std::endl;

    // auto [val63, err63] = GaussKronrodIntegral<63>(fresnel, interval);
    // std::cerr<<"GK63     = "<<std::setprecision(16)<<val63<<" "<<err63<<std::endl;

        
  }

}//namespace
