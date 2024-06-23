namespace ffd::fft{

  template<std::size_t n, direction d, template<typename, std::size_t> typename array_t>
  
  array_t<Complex, (1ul<<n)>

  FFT_l(array_t<Complex, (1ul<<n)> const& f,
	array_t<Complex, (1ul<<n)> const& phase){
    static_assert( d == direct || d == inverse );
	std::size_t constexpr pow_2_n = (1ul<<n);

    
	assert(( pow_2_n == size(f) ));

    
	array_t<Complex, pow_2_n> v_old;
	array_t<Complex, pow_2_n> v_new;
	for(std::size_t j=0; j<pow_2_n; ++j){
	  v_old[j] = f[ invert_binary_number(j, n) ];
	}

    
	for(std::size_t u = 0; u < n; ++u){
	  for(std::size_t j = 0; j < pow_2_n; ++j){
	    std::size_t const two_j = (2*j)%pow_2_n;
	    std::size_t const shift = n - u - 1;
	    std::size_t phase_index = j >> shift;
	    phase_index <<= shift;
	    v_new[j] = v_old[two_j] +
	      v_old[two_j+1]*phase[ phase_index ];
	  }
	  v_old = v_new;
	}
    
    
    
	return v_old;
      }
    
    
    
	template<std::size_t n, direction d, template<typename, std::size_t> typename array_t>
  
	array_t<Complex, (1ul<<n)>

	FFT_l(array_t<Complex, (1ul<<n)> const& f){
	  return FFT_l<n, d>(f, create_phases<(1ul<<n), d, array_t>());
	}
  template<std::size_t n, direction d>
  
  std::array<Complex, (1ul<<n)>

  FFT(std::array<Complex, (1ul<<n)> const& f,
      std::array<Complex, (1ul<<n)> const& phase){
    return FFT_l<n, d, std::array>(f, phase);
  }
    
    
    
    template<std::size_t n, direction d>
  
    std::array<Complex, (1ul<<n)>

    FFT(std::array<Complex, (1ul<<n)> const& f){
      return FFT<n, d>(f, create_phases<(1ul<<n), d>());
    }


}//namespace
