namespace ffd::chebyshev_polynomial::unit_test{

  void purify(){

    Real precision = 1e-12;

    auto f = [](auto x){ return std::sin(x*sin(x))*std::exp(-x/(1+x*x));};

    auto F = ChebyshevPolynomial<Real>(f, {-1., 2.}, precision);
    auto F_pure = F;
    // std::cerr<<std::setprecision(15)<<F(0.2)<<" "<<size(F)<<std::endl;


    F_pure.Purify(precision);

    
    int const N_samples = 100;
    for(int j=0; j<N_samples; ++j){
      Real x = -1 + 3.*(j+.2)/N_samples;
      // std::cerr<<std::abs(F(x) - F_pure(x))<<std::endl;
      assert((
	      std::abs(F(x) - F_pure(x)) <
	      2*precision
	    ));
    }

    // std::cerr<<std::setprecision(15)<<F(0.2)<<" "<<size(F)<<std::endl;



    // Real precision2 = 1e-10;
    
    // F.Purify(precision2);

    // std::cerr<<std::setprecision(15)<<F(0.2)<<" "<<size(F)<<std::endl;



    // Real precision3 = 1e-8;
    
    // F.Purify(precision3);

    // std::cerr<<std::setprecision(15)<<F(0.2)<<" "<<size(F)<<std::endl;

    
  }

}//namespace
