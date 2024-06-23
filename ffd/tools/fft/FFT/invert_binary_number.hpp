namespace ffd::fft{

  template<typename IntType>
  
  IntType
  
  invert_binary_number(IntType x, int n_digits){
    assert(( (int)x  >=  0 ));

    
    IntType x_inverted = 0;


    for(unsigned long j = 1ul<<(n_digits-1);
	x>0;
	j >>= 1, x >>= 1){
      x_inverted += j*(x%2);
    }

    
    return x_inverted;
  }


}//namespace
