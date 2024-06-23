namespace ffd::fft::unit_test{

  template<std::size_t nx, std::size_t ny>
  bool
  idempotency_2d(std::array<Complex, (1<<(nx+ny))> const& f){
    bool IsOk = true;
    std::size_t constexpr size_x_y = 1<<(nx+ny);

    
    auto fft_f = FFT2D<nx, ny>(f);
    auto fft_i = FFT2D<nx, ny, inverse>(fft_f);
    

    
    for(std::size_t j=0; j<size_x_y; ++j){
      IsOk = IsOk &&
      	std::abs( fft_i[j]/Real(size_x_y) - f[j] ) <
      	200*std::numeric_limits<Real>::epsilon();
      // std::cerr<<j<<" "<<Complex(f[j])<<" "<<fft_i[j]/Real(size_x_y)<<std::endl;
    }
    
    
    return IsOk;
  }


}//namespace
