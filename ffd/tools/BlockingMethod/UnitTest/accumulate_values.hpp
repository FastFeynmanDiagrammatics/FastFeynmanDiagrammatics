namespace ffd::blocking_method::unit_test{

  void accumulate_values(){
    using namespace ffd::user_space;

    
    block_t acc_values;

    
    for(int j: Range(297) ){
      Real is_even = ((j%2) == 0) ? 1 : -1;
      acc_values += is_even;
      // acc_values.add(is_even);
    }


    [[maybe_unused]] auto [mean, error] = acc_values.mean_error();
    // std::cerr<<mean<<" "<<error<<std::endl;
    
  }

}//namespace
