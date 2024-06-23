namespace ffd::sunset_propagator{
	
  using cheby_t = ffd::chebyshev_polynomial::ChebyshevPolynomial<Real>;

  template<ffd::phys::diagrammatic_scheme scheme,
	   int Lx, int Ly, int n_iterations = 300>
  
  std::pair<std::array<cheby_t, Lx*Ly>, std::array<cheby_t, Lx*Ly>>

  SunsetPropagator(std::array<cheby_t, Lx*Ly> const& G0,
		   Real const Beta,
		   Real const U,
		   Real const alpha_iter,
		   Real const precision = 1e-12);


  template<std::size_t n, ffd::phys::statistics stat>
  using cheby_s = ffd::chebyshev_polynomial_s::
    TPoly<n, Real, stat, ffd::l_array>;
  
  template<std::size_t n>
  using fermi_cheby_s = cheby_s<n, ffd::phys::fermi>;
  
  template<std::size_t n>
  using bose_cheby_s = cheby_s<n, ffd::phys::bose>;
  
  template<std::size_t l_x, std::size_t l_y, std::size_t n_t, ffd::phys::statistics stat>
  using green_function_t = ffd::l_array<cheby_s<n_t, stat>, (1ul<<(l_x+l_y))>;

  template<std::size_t l_x, std::size_t l_y, std::size_t n_t>
  using fermi_G_t = green_function_t<l_x, l_y, n_t, ffd::phys::fermi>;


  template<std::size_t l_x, std::size_t l_y, std::size_t n_t>
  using bose_G_t = green_function_t<l_x, l_y, n_t, ffd::phys::bose>;


  
  template<std::size_t l_x, std::size_t l_y, std::size_t n_t>
  
  fermi_G_t<l_x, l_y, n_t>

  iteration_dyson_equation(fermi_G_t<l_x, l_y, n_t> const& G0,
			   fermi_G_t<l_x, l_y, n_t> const& Sigma,
			   fermi_G_t<l_x, l_y, n_t> const& G,
			   Real const alpha_iter = 0.);


  
  template<std::size_t l_x, std::size_t l_y, std::size_t n_t>
  
  fermi_G_t<l_x, l_y, n_t>

  sunset_self_energy(fermi_G_t<l_x, l_y, n_t> const& G,
		     Real const U);
  

  
  template<std::size_t l_x, std::size_t l_y, std::size_t n_t,
	   std::size_t n_iter = 100,
	   ffd::phys::diagrammatic_scheme scheme = ffd::phys::bold>
  
  std::array<fermi_G_t<l_x, l_y, n_t>, 2>

  sunset_propagator_self_energy(fermi_G_t<l_x, l_y, n_t> const& G0,
				Real const U,
				Real const alpha_iter);
  
    
}//namespace
