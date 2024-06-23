namespace ffd::user_space::qfs_ostream::unit_test{

  void crazy_action(){
    auto den = (Bar(Psi())*Psi())('u');
    auto den_ok = den-den;
    auto den_crazy = HermitianConjugate(den*den(0.1)*Complex(1, 1))- den(-1.2);
    auto super_crazy = den_crazy*FlipSpin(den_crazy);
    std::cerr << super_crazy << '\n';
  }

}//namespace
