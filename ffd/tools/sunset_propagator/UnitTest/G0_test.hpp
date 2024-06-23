namespace ffd::sunset_propagator::unit_test{

  template<std::size_t l_x, std::size_t l_y,
	   std::size_t n_t>
  Real G0_test(Real const beta,
	       Real const mu0,
	       Real const t = 1.,
	       Real const tp = 0.){
    auto G0_kt =
      G0_k_time<l_x, l_y, n_t, ffd::l_array>(beta, mu0, t, tp);
    
    
    auto G0_rt = from_ktime_to_spacetime<l_x, l_y, n_t>(G0_kt);
    

    std::size_t const N_tau = 10;
    Real diff_max = 0.;
    for(std::size_t r=0; r<(1ul<<(l_x+l_y)); ++r){
      for(std::size_t k=0; k<N_tau; ++k){
	Real const tau = (k+.2)*beta/N_tau;
	auto g0 = -2.*G0_rt[r](tau);
	auto [x, y] = ffd::math_tools::ModDiv(r, 1ul<<l_x);
	auto g0_ex = -2*G0_help::G0_square_lattice<(1<<l_x), (1<<l_y)>(x, y, tau, beta, mu0, tp);
	Real const diff = std::abs(g0-g0_ex);
	if(diff>diff_max){
	  diff_max = diff;
	}
	// std::cerr<<g0<<" "<<diff<<" "<<std::endl;
      }
    }

    
    auto G0_kt2 = from_spacetime_to_ktime<l_x, l_y, n_t>(G0_rt);

    
    for(std::size_t K=0; K<(1ul<<(l_x+l_y)); ++K){
      for(std::size_t k=0; k<N_tau; ++k){
	Real const tau = (k+.2)*beta/N_tau;
	auto g0 = -G0_kt2[K](tau);
	auto g0_ex = -G0_kt[K](tau);
	Real const diff = std::abs(g0-g0_ex);
	if(diff>diff_max){
	  diff_max = diff;
	}
	// std::cerr<<g0<<" "<<diff<<" "<<std::endl;
      }
    }

    return diff_max;
  }

}//namespace
