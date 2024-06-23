namespace ffd::user_space::quantum_field_sum::unit_test{

  void creation(){
    auto psi = Bar(Psi('u'));
    auto psi2 = Psi('d');

    // std::cerr << std::boolalpha <<

    //   (psi[0][0].first.direction == phys::ou)
    // 	      << "\n";
    // std::cerr << psi << '\n';;
    // std::cerr << Psi() << '\n';
    // std::cerr << Psi('u') << '\n';
    // std::cerr << Bar(Psi('u')) << '\n';
  }

}//namespace
