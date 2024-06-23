namespace ffd::diagmc_sampling::unit_test{

  void creation(){
    using ffd::vector_range::Range;
    

    int const order_max = 2;
    std::vector<Real> poly;
    for( auto j: Range(1<<order_max) ){
      poly.push_back( ffd::user_space::Proba() );
    }

    
    std::array<Real, 2> seed_values{1., 2.};


    [[maybe_unused]] auto [vertex_number,
			   time_low_order,
			   time_high_order] =
      HeatBath<2, 1>(poly,
		  seed_values);


    // std::cerr<<vertex_number<<" "<<time_highest_order<<std::endl;
    

  }

}//namespace
