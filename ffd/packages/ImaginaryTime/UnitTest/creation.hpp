

namespace ffd::imaginary_time::unit_test{

  void creation(){
    
    PeriodicImaginaryTime IT(1.);

    
    for(int j=0; j<2; ++j){
      assert( std::abs(IT.LowerUpperBound[j]  - (j==1)*1.) <
	      10*std::numeric_limits<Real>::epsilon() );
    }
    
    AntiperiodicImaginaryTime ITA(2.);
    
  }

}//namespace
