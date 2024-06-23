namespace ffd::sunset_propagator{

  template<int Lx, int Ly>

  auto
  compute_sigma_sunset(std::array<cheby_t, Lx*Ly> const& G0,
		       Real const Beta,
		       Real const U,
		       Real const precision){
    std::array<cheby_t, Lx*Ly> Sigma;

    Real const U2 = U*U;
    for(int x=0; x<Lx; ++x){
      for(int y=0; y<Ly; ++y){
	int const sp_ind = y+x*Ly;
	Sigma[sp_ind] = cheby_t([=](Real tau){
				  Real const G0_tau = G0[sp_ind](tau);
				  return U2*G0_tau*G0[sp_ind](Beta-tau);},
	  {0., Beta}, precision);
	Sigma[sp_ind].Purify(precision);
      }
    }
    
    
    return Sigma;
  }

}//namespace
