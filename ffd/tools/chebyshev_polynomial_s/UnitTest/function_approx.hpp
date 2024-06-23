namespace ffd::chebyshev_polynomial_s::unit_test{

  template<std::size_t n, typename function_t>
  Real function_approx(function_t const& f){

    
    TPoly<n> G(f, {0., 1.});
    //    std::cerr << ToString(G) << '\n';
    std::stringstream ss;
    ss << ToString(G, 13) << '\n' << ToString(G, 13);
    TPoly<n> F, H;
    H.FromStream(ss);
    F.FromStream(ss);
    //    std::cerr << ToString(F) << '\n';
    std::size_t const N_samples = 1000;
    Real diff_max = 0.;
    for(std::size_t j=0; j<N_samples; ++j){
      Real const x = (j+.3)/N_samples;
      // std::cerr<<F(x)<<" "<<f(x)<<std::endl;
      Real const diff = std::abs(F(x)-f(x));
      if(diff_max<diff){
	diff_max = diff;
      }
    }
    //    std::cerr << ToString(F) << '\n';

    return diff_max;
  }

}//namespace
