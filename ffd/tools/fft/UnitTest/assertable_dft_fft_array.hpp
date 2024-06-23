namespace ffd::fft::unit_test{

  template<std::size_t n>
  bool
  assertable_dft_fft_array(std::array<Complex, (1<<n)> const& f){
    bool IsOk = true;
    using ffd::vector_range::Range;
    std::vector<Complex> f_vec(1<<n);
    for(std::size_t j=0; j<(1<<n); ++j){
      f_vec[j] = f[j];
    }
    
    
    auto dft_f = DFT(f_vec);
    auto fft_f = FFT<n>(f);
    auto fft_f_vec = FFT(f_vec);

    
    assert(( dft_f.size() == fft_f.size() ));

    
    for(  auto j: Range( fft_f.size() )  ){
      IsOk = IsOk &&
	std::abs( fft_f[j] - dft_f[j] ) <
	200*std::numeric_limits<Real>::epsilon();
      // std::cerr<<j<<" "<<fft_f[j]<<" "<<dft_f[j]<<" "<<fft_f_vec[j]<<std::endl;
    }
    
    
    return IsOk;
  }


}//namespace
