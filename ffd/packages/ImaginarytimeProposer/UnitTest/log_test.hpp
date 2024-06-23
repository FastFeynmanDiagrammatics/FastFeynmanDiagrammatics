namespace ffd::itime_proposer::unit_test{

  void
  log_test(){
    Real const Beta = 1.;
    Real const t0 = 1.;
    auto f = log_g(Beta, t0);
    std::cerr << f() << ' '
	      << f(.3) << ' '
	      << f(.2, .3) << '\n';
  }

}//namespace
