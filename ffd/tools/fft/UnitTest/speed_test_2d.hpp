namespace ffd::fft::unit_test{

  template<std::size_t nx, std::size_t ny = nx>
  void speed_test(){
    using ffd::user_space::Proba;
    std::size_t constexpr size_x = 1<<nx;
    std::size_t constexpr size_x_y = 1<<(nx+ny);

    
    std::vector<Complex> f_vec(size_x_y);
    std::array<Complex, size_x_y> f_arr;
    ffd::user_space::Timer time_dft, time_fft;

    
    for(std::size_t j=0; j<size_x_y; ++j){
      Complex const z = Complex(2.*Proba()-1., 2.*Proba()-1.);
      f_arr[j] = z;
      f_vec[j] = z;
    }
    

    time_fft.ini();
    auto const f_fft = FFT2D<nx, ny>(f_arr);
    time_fft.fin();

    
    time_dft.ini();
    auto const f_dft = DFT2D(f_vec, size_x);
    time_dft.fin();




    bool IsOk = true;
    Real max_diff = 0.;
    for(std::size_t j=0; j<size_x_y; ++j){
      Real diff_now = std::abs(f_fft[j] - f_dft[j]);
      if(diff_now > max_diff){
	max_diff = diff_now;
      }
      IsOk = IsOk &&
	diff_now < 500*std::numeric_limits<Real>::epsilon();
    }
    // if(IsOk){
    //   std::cerr<<"DFT2D vs FFT2D: same result obtained\n";
    // }else{
    //   std::cerr<<"DFT2D vs FFT2D: different result obtained !!!!!!!!!!!!!!!!!!!!!!!!!\n";
    // }
    std::cerr<<"max_difference(DFT2D-FFT2D) = "<<max_diff<<"\n";
    

    std::cerr<<"time(DFT2D) = "<<time_dft()<<"s, time(DFT2D)/time(FFT2D) = "<<time_dft()/time_fft()<<"\n";
  }

}//namespace
