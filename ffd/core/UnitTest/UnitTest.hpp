#ifdef FFD_UNIT_TEST_FLAG

namespace ffd{

  int counter_unit_test = 0;

  template<typename test_t>
  void unit_test_runner(test_t test_v){
    test_v();

    
    std::ios_base::fmtflags flags_cerr_saved( std::cout.flags() );  // save flags state
    std::cout<<std::hex<<counter_unit_test<<" "<<std::flush;
    std::cout.flags( flags_cerr_saved );
    ++counter_unit_test;
  }

}//namespace

namespace ffd::core_unit_test{

  void UnitTest(){
    
    
    ffd::counter_unit_test = 0;

    std::cerr<<"test# = ";
    
    unit_test_runner(ffd::core_math::unit_test::UnitTest);
    
    unit_test_runner(ffd::set_theory::unit_test::UnitTest);

    unit_test_runner(ffd::flat_map::unit_test::UnitTest);

    unit_test_runner(ffd::nilpotent_polynomial::unit_test::UnitTest);
    
    unit_test_runner(ffd::ranked_zeta::unit_test::UnitTest);

    unit_test_runner(ffd::quantum_field::unit_test::UnitTest);
  
    unit_test_runner(ffd::qf_dot::unit_test::UnitTest);
    
    unit_test_runner(ffd::qf_product::unit_test::UnitTest);

    unit_test_runner(ffd::qf_vertex::unit_test::UnitTest);

    unit_test_runner(ffd::qf_graph::unit_test::UnitTest);

    unit_test_runner(ffd::feynman_edge::unit_test::UnitTest);

    unit_test_runner(ffd::wick_function::unit_test::UnitTest);

    unit_test_runner(ffd::wick_matrix::unit_test::UnitTest);

    std::cerr<<" ||  ";


  }

}

#endif
