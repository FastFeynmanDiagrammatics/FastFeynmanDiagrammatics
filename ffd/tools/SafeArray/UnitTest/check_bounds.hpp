namespace ffd::user_space::safe_array::unit_test{

  void check_bounds(){
    
    SafeArray<Real, 2> x;

    ffd::get<0>(x) = 3.14;

    ffd::get<1>(x) = 6.28;

    // ffd::get<2>(x) = 6.28;

  }


}//namespace
