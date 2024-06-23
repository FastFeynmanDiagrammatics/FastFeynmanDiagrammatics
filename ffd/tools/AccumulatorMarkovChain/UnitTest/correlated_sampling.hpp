namespace ffd::accumulator_markov_chain::unit_test{

  void correlated_sampling(){
    using ffd::random_distributions::Proba;
    AccumulatorMarkovChain G(100);
    ffd::binary_accumulator::BinaryAccumulator<Real> B;

    long const iterations = 1000ul;
    
    Real x_now = 1;
    Real const x_mean = -.2;
    Real const alpha = 0.99;
    Real x_acc = 0;
    for(int j=0; j < iterations; ++j){
      x_now = alpha*(x_now-x_mean) + 2*(Proba()-.5) + x_mean;
      G += x_now;
      B += x_now;
      x_acc += x_now;
    }
    
    // std::cerr<<G.Variance()<<std::endl;
    // for(int j=1; j<100; ++j){
    //   std::cerr<<G.AutocorrelationFunction(j)<<std::endl;
    // }
    // std::cerr<<G.AutocorrelationTime()<<" "<<alpha/(1-alpha)<<std::endl;

    [[maybe_unused]] auto [value, error] = G.MeanValueError();

    // std::cerr<<x_acc/iterations<<std::endl;
    // std::cerr<<B.MeanValue()<<std::endl;
    // std::cerr<<value<<" +/- "<<error<<std::endl;
    // std::cerr<<x_mean<<" dif "<<std::abs(value-x_mean)<<std::endl;
  }

}//namespace
