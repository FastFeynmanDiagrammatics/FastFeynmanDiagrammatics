namespace ffd::lattice_proposer::unit_test{

  void proposer(){
    constexpr int d = 2;
    std::array<int, d> L{{4, 5}};
    std::array<std::array<Real, d>, d> bravais;
    bravais[0] = {1., 0.};
    bravais[1] = {0., 1.};
    auto Realspace =
      ffd::lattice_g::
      Realspace_g<2>(L,
		     bravais);
      auto f =
	g<2>({4, 4},
	   [](auto x){return 1./(1.+x*x*x*x);},
	     ffd::lattice_g::Distance_g<d>(Realspace));

      ffd::lattice_g::coord_d<d> x0, x1;
      x0 = {1, 0};
      x1 = {0, 0};
      std::cerr << f(x0, x1) << '\n';
      auto x2 = f(x0);
      auto x3 = f();
      std::cerr << f(x2, x3) << '\n';

  }

}//namespace
