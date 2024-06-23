namespace ffd::fft::unit_test{

  void fft2d_l(){

    ffd::l_array<Complex, 4> v(Complex(2., 1.));


    auto x = FFT2D_l<1, 1>(v);
    

  }

}//namespace
