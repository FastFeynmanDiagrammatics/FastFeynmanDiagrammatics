namespace ffd::sunset_propagator::unit_test{


  template<int Lx, int Ly = Lx>
  [[maybe_unused]] auto
  bare_sunset(Real const Beta,
		   Real const mu,
		   Real const U,
		   Real const tp=0.,
		   Real const absolute_precision = 1e-12){
    Real const mu0 = compute_at_fixed_chemical_potential<Lx, Ly>(Beta, mu, U, tp); 
    auto const G0 = G0_cheby_square_lattice<Lx, Ly>(Beta, mu0, tp, absolute_precision);
    std::cerr<<"G0 built"<<std::endl;
    auto const [G, Sigma] = SunsetPropagator<ffd::phys::bare, Lx, Ly, 1>(G0, Beta, .8, absolute_precision);



    int const N_tau = 1ul<<10;
    Real const Pi = 3.1415926535;
    Real const omega = Pi/Beta;
    for(int kx=0; kx<=Lx/2; ++kx){
      Real const Kx = kx*2.*Pi/Lx;
      for(int ky=0; ky<=Lx/2; ++ky){
	Real const Ky = ky*2.*Pi/Ly;
	Real im_s = 0., re_s = 0.;
	for(int x=0; x<Lx; ++x){
	  for(int y=0; y<Ly; ++y){
	    for(int j=0; j<N_tau; ++j){
	      Real const tau = (.5+j)*Beta/N_tau;
	      Real const common_part = cos(Kx*x)*cos(Ky*y)*Sigma[y+x*Ly](tau);
	      im_s += sin(omega*tau)*common_part;
	      re_s += cos(omega*tau)*common_part;
	    }
	  }
	}
	im_s /= Lx*Ly*N_tau;
	re_s /= Lx*Ly*N_tau;
	std::cerr<<kx<<" "<<ky<<" "<<re_s<<" "<<im_s<<"\n";
      }
    }

    
    std::cerr<<-2.*G[0](Beta)<<std::endl;

    
    return std::make_pair(G, Sigma);
  }

}//namespace
