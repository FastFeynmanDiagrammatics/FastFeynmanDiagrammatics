namespace ffd::fft::unit_test{

  template<std::size_t n>
  bool
  idempotency(std::array<Complex, (1<<n)> const& f){
    bool IsOk = true;
    std::size_t constexpr two_n = 1<<n;
    
    auto fft_f = FFT<n>(f);
    auto fft_i = FFT<n, inverse>(fft_f);

    

    
    for(std::size_t j=0; j<two_n; ++j){
      IsOk = IsOk &&
      	std::abs( fft_i[j]/Real(two_n) - f[j] ) <
      	200*std::numeric_limits<Real>::epsilon();
      // std::cerr<<j<<" "<<Complex(f[j])<<" "<<fft_i[j]/Real(two_n)<<std::endl;
    }
    
    
    return IsOk;
  }




}//namespace
