namespace ffd::fft::unit_test{

  void fft_l(){
    
    ffd::l_array<Complex, 2> v(1.);
    
    auto f = FFT_l<1>(v);
    
  }
  
}//namespace
