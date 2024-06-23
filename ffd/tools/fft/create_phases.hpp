namespace ffd::fft{

  auto
  create_phases(unsigned long pow_2_n){
    std::vector<Complex> phases(pow_2_n);

    
    Real two_pi_over_2_n = 2*ffd::core_math::Pi/pow_2_n;
    for( unsigned long j=0; j<pow_2_n; ++j ){
      phases[j] = ffd::user_space::expI( j*two_pi_over_2_n*(-1) );
    }

    
    return phases;
  }

  

  template<std::size_t pow_2_n,
	   direction d = direct,
	   template<typename, std::size_t> typename array_t = std::array>
  
  constexpr
  array_t<Complex, pow_2_n>
  create_phases(){
    static_assert( d == direct || d == inverse );
    array_t<Complex, pow_2_n> phases;


    Real constexpr two_pi_over_2_n = 2.*ffd::core_math::Pi*(-1.+2*(d==inverse))/pow_2_n;
    for( unsigned long j=0; j<pow_2_n; ++j ){
      phases[j] = ffd::user_space::expI( j*two_pi_over_2_n );
    }

    
    return phases;
  }


}//namespace
