

namespace ffd::contour_derivative::unit_test{

  void one_over_one_minus_x(){

    
    auto f = [](auto x){
	       return 1./(1.-x);
	     };


    ContourDerivative fd(f, .3);

    using ffd::vector_range::Range;
    for( auto j: Range(12) ){
      // std::cerr<<fd(j)<<std::endl;

      assert((
	      std::abs( fd(j)/ffd::core_math::Factorial(j) - 1.) < 1e-10
	      ));
    }

    
    
  }


}//namespace
