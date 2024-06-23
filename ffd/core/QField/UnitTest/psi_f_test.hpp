namespace ffd::user_space::q_field::unit_test{

  void
  psi_f_test(){
    auto p = psi_f(1);
    std::cerr << p << '\n';
    std::cerr << Component(p) << '\n';
    p = psi_f(-1);
    std::cerr << p << '\n';
    std::cerr << Component(p) << '\n';
    std::cerr << sizeof(p) << '\n';
  }

}//namespace
