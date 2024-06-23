namespace ffd::itime_lattice_proposer::unit_test{

  void itime_l(){
    auto f =
      [](auto ...x){return 2;};
    auto h =
      [](auto ...x){return -1;};
	
    auto fh = g(f, h);

    auto [x, y] = fh(std::make_pair(1, 2));
    std::cerr << x << ' ' << y << '\n';
  }

}//namespace
