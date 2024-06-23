namespace ffd::sigmadet_hubbard::unit_test{

  void UnitTest(){

    instantiation<2>();

    instantiation<3>();

    // instantiation<4>();
    CDet_old_vs_new<2>();

    CDet_old_vs_new<7>();
    CDet_old_vs_new<8>();
    
  }

}//namespace
