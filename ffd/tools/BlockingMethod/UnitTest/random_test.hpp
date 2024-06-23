namespace ffd::blocking_method::unit_test{

  void random_test(){
    using namespace ffd::user_space;
    
    
    block_t accumulator;
    Real const alpha = 0.7;
    Real random_old = 0.;

    
    for( int j: Range( 121 ) ){
      Real random_real = 2*Proba() - 1.;
      Real random_corr = random_old*alpha + (1-alpha)*random_real;
      random_old = random_corr;
      // accumulator.add(random_real);
      accumulator += random_corr;
    }
    

    // [[maybe_unused]] auto [mean, error] = accumulator.mean_error();
    [[maybe_unused]] auto [mean, error, auto_time, variance] = accumulator.statistics();

    // std::cerr<<mean<<" "<<error<<" "<<auto_time<<" "<<variance<<std::endl;
    // std::cerr<<"0 nan "<<alpha/(1-alpha)<<"\n";
    
  }


}//namespace
