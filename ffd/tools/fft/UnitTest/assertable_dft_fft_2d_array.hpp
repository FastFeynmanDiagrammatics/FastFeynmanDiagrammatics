namespace ffd::fft::unit_test{

  template<std::size_t nx,
	   std::size_t ny>
  bool
  assertable_dft_fft_2d_array(std::array<Complex, (1<<(nx+ny))> f){
    bool IsOk = true;
    using ffd::vector_range::Range;
    std::size_t constexpr size_x = (1<<nx);
    std::size_t constexpr size_y = (1<<ny);
    std::size_t constexpr size_x_y = (1<<(nx+ny));
    std::vector<Complex> f_vec(size_x_y);
    for(std::size_t r=0; r<size_x_y; ++r){
      f_vec[r] = f[r];
    }
    
    
    auto dft_f = DFT2D(f_vec, size_x);
    auto fft_f_vec = FFT2D(f_vec, size_x);
    auto fft_f = FFT2D<nx, ny>(f);

    
    assert(( dft_f.size() == fft_f.size() ));

    
    for(  auto j: Range( size(fft_f) )  ){
      IsOk = IsOk &&
	std::abs( fft_f[j] - dft_f[j] ) <
	200*std::numeric_limits<Real>::epsilon();
      // std::cerr<<j<<" "<<fft_f[j]<<" "<<dft_f[j]<<" "<<fft_f_vec[j]<<std::endl;
    }
    
    
    return IsOk;
  }


}//namespace
