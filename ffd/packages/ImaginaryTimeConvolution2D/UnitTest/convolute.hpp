namespace ffd::imaginary_time_convolution_2d::unit_test{

  void convolute(){
    Real a = 2.12313, b=1.12312, c=0.8213124;
    Real beta = 2.12312;
    Real precision, precision_fast;
    
    auto f = [a](Real){return a;};
    auto g = [b](Real){return b;};
    auto h = [c](Real){return c;};
    bool g_is_fermion, h_is_fermion;


    auto L = Convolute123FromZeroTo4<Real>(f, g, h,
					   beta,
					   {{g_is_fermion=true,
						 h_is_fermion=true}},
					   precision=.02);



    //it provides the right answer for tau1<=tau2,
    //otherwise it gives an analytic continuation
    auto L_fast = FastFermionicConvolute123FromZeroTo4<Real>(f, g, h,
							     beta,
							     precision_fast=1e-10);

    auto L_very_fast = VeryFastFermionicConvolute123FromZeroTo4Dflt(f, g, h,
								    beta,
								    precision_fast=1e-10);

    


    auto L_exact = [beta, a, b, c](Real tau1, Real tau2){
		     return a*b*c*(beta - 2*std::abs(tau1-tau2));
		   };

    

    // std::cerr<<"order = "<<L.Order<<std::endl;

    // std::cerr<<"order_fast = "<<L_fast.Order<<std::endl;

    const int N_samples = 10;
    for(int y=0; y < N_samples; ++y){
      for(int x=0; x <= y; ++x){
	Real xx = (x+.1231231)*beta/N_samples;
	Real yy = (y+.1231231)*beta/N_samples;
	// std::cerr<<"("<<xx<<", "<<yy<<") "<<(L_fast(xx, yy)-L_exact(xx, yy))<<
	//   " "<<(L(xx, yy)-L_exact(xx, yy))<<" VF "<<(L_very_fast(xx, yy)-L_exact(xx, yy))<<std::endl;
	assert( std::abs(L(xx, yy)-L_exact(xx, yy)) < precision);
	assert( std::abs(L_fast(xx, yy)-L_exact(xx, yy)) < precision_fast);
	assert( std::abs(L_very_fast(xx, yy)-L_exact(xx, yy)) < precision_fast);
      }
    }
    
  }


}//namespace
