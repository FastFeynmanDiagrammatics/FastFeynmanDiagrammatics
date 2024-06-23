namespace ffd::sunset_propagator::unit_test{

  template<int Lx, int Ly = Lx>
  [[maybe_unused]]
  Real
  compute_at_fixed_chemical_potential(Real const Beta,
				      Real const mu,
				      Real const U,
				      Real const tp=0.){

    auto mu0 = ffd::find_root::FindRealRoot([=](Real mu_0){return
	  mu_0 - G0_help::G0_square_lattice<Lx, Ly>(0, 0, Beta, Beta, mu_0, tp)*U-mu;},
      {-10., mu});


    std::cerr<<std::setprecision(15)<<
      "n0 = "<<-2*G0_help::G0_square_lattice<Lx, Ly>(0, 0, Beta, Beta, mu0.value(), tp)
	     <<", mu0 = "<<mu0.value()<<std::endl;

    return mu0.value();
  }

}//namespace
