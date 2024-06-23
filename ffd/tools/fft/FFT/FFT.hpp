namespace ffd::fft{

  std::vector<Complex>

  FFT(std::vector<Complex> f){
    using ulong_t = unsigned long;
    

    ulong_t const pow_2_n = f.size();
    ulong_t const log2_size = ffd::math_tools::log2_int(pow_2_n);
    assert((
	    ffd::core_math::pow_int(2, log2_size) == (long)pow_2_n
	    ));

    

    auto phase = create_phases(pow_2_n);


	    
    std::vector<Complex> v_old(pow_2_n), v_new(pow_2_n);
    for( ulong_t j=0; j<pow_2_n; ++j){
      v_old[j] = f[ invert_binary_number(j, log2_size) ];
    }


    
    for( ulong_t u = 0; u < log2_size; ++u){
      for( ulong_t j = 0; j < pow_2_n; ++j){
	ulong_t two_j = (2*j)%pow_2_n;
	int const shift = log2_size - u - 1;
	ulong_t phase_index = j >> shift;
	phase_index <<= shift;
	v_new[j] = v_old[two_j] +
	  v_old[two_j+1]*phase[ phase_index ];
      }
      v_old = v_new;
    }
    

    
    return v_old;
  }


}//namespace
