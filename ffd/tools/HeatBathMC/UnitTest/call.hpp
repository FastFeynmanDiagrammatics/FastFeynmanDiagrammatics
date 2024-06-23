namespace ffd::heat_bath_mc::unit_test{

  void UnitTest(){

    partial_sum_test();
    //std::cout << "passed PS" << std::endl;   
    conformal_map_test();  
    //std::cout << "passed CM" << std::endl;
    //exit(1);

    int constexpr order = 4;
    int constexpr order_plus = 1;
    int constexpr half_sites = 4;
    unsigned long constexpr N_iter = 1ul<<23;
    
    // exact_sites<order, order_plus, 1>();
    // exact_sites<order, order_plus, 1, false>();
    // exact_sites<order, order_plus, 2>();
    // exact_sites<order, order_plus, 2, false>();
    // exact_sites<order, order_plus, 3>();
    // exact_sites<order, order_plus, 3, false>();

    
    invoke_exact_all_orders<order, order_plus, half_sites>(10., N_iter);
    invoke_translation_invariant<order, order_plus, half_sites>(.33, N_iter, 9.1231);
  }

}//namespace
