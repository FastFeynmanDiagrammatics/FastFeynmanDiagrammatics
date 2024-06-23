namespace ffd::fft{

  std::vector<Complex>
  FFT2D(std::vector<Complex> f,
	unsigned long size_x = 0ul);


  
  template<std::size_t nx,
	   std::size_t ny = nx,
	   direction d = direct>

  std::array<Complex, (1ul<<(nx+ny))>

  FFT2D(std::array<Complex, (1ul<<(nx+ny))> const&);

  

  // template<std::size_t nx,
  // 	   std::size_t ny = nx,
  // 	   direction d = direct,
  // 	   template<typename, std::size_t> typename array_t = std::array>

  // array_t<Complex, (1ul<<(nx+ny))>

  // FFT2D_l(array_t<Complex, (1ul<<(nx+ny))> const&);


}//namespace
