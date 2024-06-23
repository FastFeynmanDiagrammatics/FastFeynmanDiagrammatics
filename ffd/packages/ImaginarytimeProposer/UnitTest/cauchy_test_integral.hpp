namespace ffd::itime_proposer::unit_test{

  void
  cauchy_test_integral(){
    Real const t0 = 1.2122131, t1 = 2.12312312, t2 = 1.43234234;
    auto f = cauchy_g(t0);
    auto g = [t1](Real t){return .5*exp(-std::abs(t/t1))/t1;};
    long MC_iter = 1ul<<28;
    Real avg = 0;
    for(std::size_t j=0; j<MC_iter; ++j){
      auto t = f(t2);
      avg += g(t)/f(t2, t);
    }
    std::cerr << "0 = " << std::setprecision(10) << avg/MC_iter-1. << '\n';
  }

}//namespace
