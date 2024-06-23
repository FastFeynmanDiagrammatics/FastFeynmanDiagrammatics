namespace ffd::fft::unit_test{

  void UnitTest(){

    check_invert_binary_number();
    
    compare_dft_fft();

    compare_dft_fft_2d();

    compare_fft_array_vector();

    check_idempotency();

    compare_dft_fft_2d_array();

    // speed_test<8>();

    check_idempotency_2d();

    fft_l();
    
  }

}//namespace
