namespace ffd::chebyshev_polynomial_s::unit_test{

  template<std::size_t n, typename field_t = Real>
  Real convolution(Real const Beta){
    
    auto const f = [](Real x){return
	// 1./(1+x+x*x);};
	std::sin(3*x*x)/(1+x*x);};
    auto const g = [](Real x){return
	// 1./(1+x*x+x*x*x*x);};
	std::cos(5*x*x/(1+x*x))/(1.+x);};
    auto const Tf = TPoly<n, field_t>(f, {0., Beta});
    auto const Tg = TPoly<n, field_t>(g, {0., Beta});
    
    
    auto const Th = convolute(Tf, Tg);


    int const N_samples = 200;
    Real diff_max = 0.;
    for(int j=0; j<N_samples; ++j){
      Real const t = (j+.2)*Beta/N_samples;
      Real const conv1 = ffd::gauss_kronrod::
	Integrate_63([=](Real t1){return f(t1)*g(t-t1);},
		     {0., t});
      Real const conv2 = ffd::gauss_kronrod::
	Integrate_63([=](Real t1){return -f(t1)*g(Beta+t-t1);},
		     {t, Beta});
      auto const Th_t = Th(t);
      Real const diff = std::abs(conv1+conv2-Th(t));
      if(diff>diff_max){
	diff_max = diff;
      }
      // std::cerr<<Th_t<<" "<<diff<<std::endl;
    }

    
    return diff_max;
  }

}//namespace
