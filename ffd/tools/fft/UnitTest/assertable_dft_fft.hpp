namespace ffd::fft::unit_test{

  bool
  assertable_dft_fft(std::vector<Complex> f){
    bool IsOk = true;
    using ffd::vector_range::Range;

    
    auto dft_f = DFT(f);
    auto fft_f = FFT(f);

    
    assert(( dft_f.size() == fft_f.size() ));

    
    for(  auto j: Range( fft_f.size() )  ){
      IsOk = IsOk &&
	std::abs( fft_f[j] - dft_f[j] ) <
	200*std::numeric_limits<Real>::epsilon();
      // std::cerr<<j<<" "<<fft_f[j]<<" "<<dft_f[j]<<std::endl;
    }
    
    
    return IsOk;
  }

}//namespace
