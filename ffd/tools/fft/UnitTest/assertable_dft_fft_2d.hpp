namespace ffd::fft::unit_test{

  bool
  assertable_dft_fft_2d(std::vector<Complex> f,
			unsigned long size_x = 0ul){
    bool IsOk = true;
    using ffd::vector_range::Range;

    
    auto dft_f = DFT2D(f, size_x);
    auto fft_f = FFT2D(f, size_x);

    
    assert(( dft_f.size() == fft_f.size() ));

    
    for(  auto j: Range( size(fft_f) )  ){
      IsOk = IsOk &&
	std::abs( fft_f[j] - dft_f[j] ) <
	200*std::numeric_limits<Real>::epsilon();
      // std::cerr<<j<<" "<<fft_f[j]<<" "<<dft_f[j]<<std::endl;
    }
    
    
    return IsOk;
  }



}//namespace
