namespace ffd::sunset_propagator{

  template<int Lx, int Ly>

  auto
  compute_fermionic_convolution(std::array<cheby_t, Lx*Ly> const& f,
				std::array<cheby_t, Lx*Ly> const& g,
				Real const Beta,
				Real const precision){
    std::array<cheby_t, Lx*Ly> f_g;
    
    
    for(int x=0; x<Lx; ++x){
      for(int y=0; y<Ly; ++y){
	f_g[y+x*Ly] = cheby_t([=](Real tau){
				Real val = 0.;
				for(int xx=0; xx<Lx; ++xx){
				  for(int yy=0; yy<Ly; ++yy){
				    int const sp_ind = yy+xx*Ly;
				    int const dy = (y >= yy ? y-yy : y-yy+Ly);
				    int const dx = (x >= xx ? x-xx : x-xx+Lx);
				    int const dsp_ind = dy+Ly*dx;
				    val += ffd::gauss_kronrod::
				      Integrate([=](Real t1){return f[sp_ind](t1)*g[dsp_ind](tau-t1);},
						{0., tau}, precision);
				    val -= ffd::gauss_kronrod::
				      Integrate([=](Real t1){return f[sp_ind](t1)*g[dsp_ind](Beta+tau-t1);},
						{tau, Beta}, precision);
				  }
				}
				return val;
			      }, {0., Beta}, precision);
	f_g[y+x*Ly].Purify(precision);
      }
    }


    return f_g;
  }

}//namespace
