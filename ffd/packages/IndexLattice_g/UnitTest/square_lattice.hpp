namespace ffd::index_lattice::unit_test{

  void square_lattice(){
    using namespace ffd::user_space;
    using int_d = lattice_g::int_d;
    std::array<int_d, 2> O = {0, 0}; 
    auto H0 = Bar(psi('u')(O))*psi('u')(O);
    H0 += FlipSpin(H0);
    H0 += Bar(psi(0))(O)*psi('u')(O);
    auto oracle = BlockOracle_g(H0);
    int const Lx = 3;
    int const Ly = 3;
    std::array<int, 2> L = {Lx, Ly};
    auto [indexer, tot_size]
      = g<2>(L, oracle);
    std::cerr
      << "tot_size = "
      << tot_size
      << '\n';
    {
      int_d
	x = 0,
	y = 0;
      auto X0 = std::make_pair(psi_f('d'),
			       std::array<int_d, 2>{x, y});
      auto X1 = std::make_pair(Bar(psi_f('d')),
			       std::array<int_d, 2>{0, 0});
      std::cerr
	<< X0
	<< ' '
	<< X1
	<< ' '
	<< indexer(X0, X1)
	<< '\n';
    }      
  }

}//namespace
