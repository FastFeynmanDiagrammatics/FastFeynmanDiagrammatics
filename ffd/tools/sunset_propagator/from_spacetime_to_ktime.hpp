namespace ffd::sunset_propagator{
  
  template<std::size_t l_x,
	   std::size_t l_y,
	   std::size_t n_t,
	   ffd::phys::statistics stat_zeta = ffd::phys::fermi,
	   template<typename, std::size_t> typename array_t = ffd::l_array>

  array_t<cheby_s<n_t, stat_zeta>, (1ul<<(l_x+l_y))>
  
  from_spacetime_to_ktime(array_t<cheby_s<n_t, stat_zeta>, (1ul<<(l_x+l_y))> const& G){
    std::size_t constexpr size_x_y = (1ul<<(l_x+l_y));

    
    array_t<cheby_s<n_t, stat_zeta>, size_x_y> G_kt;
    G_kt.fill(ffd::make_array(0., G[0].bounds_average*2.));
    array_t<Complex, size_x_y> G_coef;
    for(std::size_t k=0; k<n_t; ++k){
      for(std::size_t r=0; r<size_x_y; ++r){
	G_coef[r] = G[r].coef[k];
      }

      
      auto const G_coef_k = ffd::fft::FFT2D_l<l_x, l_y, ffd::fft::inverse>(G_coef);
      for(std::size_t K=0; K<size_x_y; ++K){
	G_kt[K].coef[k] = real(G_coef_k[K]);
      }
    }
    
    
    return G_kt;
  }

  

}//namespace
