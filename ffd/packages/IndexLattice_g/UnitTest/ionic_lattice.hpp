namespace ffd::index_lattice::unit_test{

  void ionic_lattice(){
    using namespace ffd::user_space;
    using int_d = lattice_g::int_d;
    ffd::lattice_g::coord_d<2, 2> O = {0, 0, 0}; 
    auto H0 = Bar(psi('u')(O))*psi('u')(O);
    H0 += FlipSpin(H0);
    //    H0 += Bar(psi(0))(O)*psi('u')(O);
    auto oracle = BlockOracle_g(H0);
    int const Lx = 2;
    int const Ly = 4;
    std::array<int, 2> L = {Lx, Ly};
    auto [indexer, tot_size]
      = g<2, 2>(L, oracle);
    std::cerr
      << "tot_size = "
      << tot_size
      << '\n';
    {
      int_d
	x0 = 1,
	y0 = 1,
	x1 = 1,
	y1 = 1,
	a0 = 1,
	a1 = 1;
      auto X0 = std::make_pair(Bar(psi_f('u')),
			       std::array<int_d, 3>{x0, y0, a0});
      auto X1 = std::make_pair(Bar(psi_f('u')),
			       std::array<int_d, 3>{x1, y1, a1});
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
