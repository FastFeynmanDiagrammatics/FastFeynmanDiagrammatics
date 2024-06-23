namespace ffd::sunset_propagator{

  template<ffd::phys::diagrammatic_scheme scheme,
	   int Lx, int Ly, int n_iterations>
  
  std::pair<std::array<cheby_t, Lx*Ly>, std::array<cheby_t, Lx*Ly>>

  SunsetPropagator(std::array<cheby_t, Lx*Ly> const& G0,
		   Real const Beta,
		   Real const U,
		   Real const alpha_iter,
		   Real const precision){
    static_assert(scheme == ffd::phys::bare || scheme == ffd::phys::bold);
    using propagator_t = std::array<cheby_t, Lx*Ly>;
    propagator_t Sigma = compute_sigma_sunset<Lx, Ly>(G0, Beta, U, precision);
    std::cerr<<"sigma_computed\n";
    propagator_t G0Sigma = compute_fermionic_convolution<Lx, Ly>(G0, Sigma, Beta, precision);
    std::cerr<<"G0Sigma computed\n";
    propagator_t G = G0;
    
    
    for(int iter=0; iter<n_iterations; ++iter){
      propagator_t G0SigmaG = compute_fermionic_convolution<Lx, Ly>(G0Sigma, G, Beta, precision);
      std::cerr<<"G0SigmaG computed\n";
      for(int r=0; r<Lx*Ly; ++r){
	G[r] += G[r]*alpha_iter + G0SigmaG[r]*(1-alpha_iter);
      }
    }


    return std::make_pair(G, Sigma);
  }


}//namespace
