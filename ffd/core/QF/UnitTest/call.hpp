namespace ffd::user_space::qf::unit_test{

  void UnitTest(){

    QF x;
    x.component = 1;
    std::cerr << x ;
    
    std::cerr << FlipSpin(x) << ' ' << Bar(x) << '\n';
  }

}//namespace
