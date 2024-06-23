namespace ffd::periodic_numbers::unit_test{

  void periodic_int(){
    auto periodize = ffd::functional::
      Ret_g<0>(Periodic_g<int>({0, 1}));
    assert((periodize(0)==0));
    assert((periodize(1)==0));
    assert((periodize(-1)==0));
    
  }

}//namespace
