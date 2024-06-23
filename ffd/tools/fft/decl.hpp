namespace ffd::fft{

  std::vector<Complex>
  FFT(std::vector<Complex> f);

  
  enum direction{direct, inverse};




  template<std::size_t n,
	   direction d = direct>
  
  std::array<Complex, (1ul<<n)>
  
  FFT(std::array<Complex, (1ul<<n)> const& f,
      std::array<Complex, (1ul<<n)> const& phases);

  
  
  template<std::size_t n,
	   direction d = direct>
  
  std::array<Complex, (1ul<<n)>
  
  FFT(std::array<Complex, (1ul<<n)> const& f);


  template<std::size_t n,
	   direction d = direct,
	   template<typename, std::size_t> typename array_t = std::array>
  
  array_t<Complex, (1ul<<n)>
  
  FFT_l(array_t<Complex, (1ul<<n)> const& f,
      array_t<Complex, (1ul<<n)> const& phases);


  template<std::size_t n,
	   direction d = direct,
	   template<typename, std::size_t> typename array_t = std::array>
  
  array_t<Complex, (1ul<<n)>
  
  FFT_l(array_t<Complex, (1ul<<n)> const& f);

	
}//namespace
