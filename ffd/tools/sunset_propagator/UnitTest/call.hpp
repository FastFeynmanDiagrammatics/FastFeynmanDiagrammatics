namespace ffd::sunset_propagator::unit_test{

  void UnitTest(){
    Real const Beta = 5.;
    Real const mu = 1.9;
    Real const tp = -.3;
    Real const U = 5.6;
    Real const precision = 1e-4;
    std::size_t constexpr l_x = 8;
    std::size_t constexpr l_y = 8;
    std::size_t constexpr n_t = 40;
    std::size_t constexpr L = 1<<l_x;
    compute_at_fixed_chemical_potential<L>(Beta, mu, U, tp);
    // auto [G, Sigma] = bare_sunset<L>(Beta, mu, U, tp, precision);

    
    // std::cerr<<(G0_test<l_x, l_y, n_t>(Beta, -.3, 1., 0.))<<std::endl;;
    assert(( G0_test<l_x, l_y, n_t>(Beta, -.5, 1., -.3) < 1e-10 ));
    
    
  }

}//namespace
