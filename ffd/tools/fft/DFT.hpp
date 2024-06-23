

namespace ffd::fft{

  std::vector<Complex>

  DFT(std::vector<Complex> f){
    unsigned long const pow_2_n = size(f);
    [[maybe_unused]] unsigned long const log2_size = ffd::math_tools::log2_int( pow_2_n );

    
    assert(( (unsigned long)ffd::core_math::pow_int(2, log2_size) == pow_2_n ));


    auto phase = create_phases(pow_2_n);


    std::vector<Complex> tilde_f(pow_2_n, 0.);
    for( unsigned long k=0; k < pow_2_n; ++k ){
      for( unsigned long j=0; j < pow_2_n; ++j ){
	auto j_k = (j*k)%pow_2_n;
	tilde_f[k] += phase[j_k]*f[j];
      }
    }

    
    return tilde_f;
  }


}//namespace
