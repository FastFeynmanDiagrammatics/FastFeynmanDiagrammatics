namespace ffd::user_space::q_f_sum::unit_test{

  void psi_test(){
    auto p = psi(1);
    auto p2 = p(std::make_pair(1, 2));
    auto p12 = Bar(p)*p*p2;
    std::cerr << p << '\n';
    std::cerr << Component(p.fields[0].first) << '\n';
    std::cerr << HermitianConjugate(Complex(1,1)*p2*FlipSpin(p2)) << '\n';
    auto S0 = .1*p12-.3*p*p12;
    std::cerr << S0 << '\n';
    std::cerr << HermitianConjugate(S0) << '\n';
  }

}//namespace
