namespace ffd::sunset_propagator::G0_help{

  inline  Real
  epsilon_k(Real kx, Real ky, Real mu0, Real tp){
    return -2*(cos(kx)+cos(ky)) - 4*tp*cos(kx)*cos(ky) - mu0;
  }
  
  
  template<int Lx, int Ly>

  auto
  G0_square_lattice(int x, int y, Real tau, Real beta, Real mu0, Real tp){
    Real ret = 0;
    
    
    using ffd::core_math::Pi;
    using std::cos, std::exp;

    
    for(int kx=0; kx<Lx; ++kx){
      for(int ky=0; ky<Ly; ++ky){
	Real k_x = 2*kx*Pi/Lx;
	Real k_y = 2*ky*Pi/Ly;
	Real Ek = epsilon_k(k_x, k_y, mu0, tp);
	ret -= cos(k_x*x)*cos(k_y*y)*exp(-tau*Ek)/(1+exp(-beta*Ek));
      }
    }
    return ret/Lx/Ly;
  }

}//namespace



namespace ffd::sunset_propagator{

  template<int Lx, int Ly>
  
  std::array<cheby_t, Lx*Ly>

  G0_cheby_square_lattice(Real const Beta,
	   Real const mu0,
	   Real const tp=0.,
	   Real const absolute_precision=1e-12){
    std::array<cheby_t, Lx*Ly> G0;


    for(int x=0; x<Lx; ++x){
      for(int y=0; y<Ly; ++y){
	G0[y+x*Ly] = cheby_t([=](Real tau){return
	      G0_help::G0_square_lattice<Lx, Ly>(x, y, tau, Beta, mu0, tp);},
	  {0., Beta}, absolute_precision);
	G0[y+x*Ly].Purify(absolute_precision);
      }
    }


    return G0;
  }

}//namespace
