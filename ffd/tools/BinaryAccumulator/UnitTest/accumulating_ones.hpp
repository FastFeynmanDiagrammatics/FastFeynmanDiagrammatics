

namespace ffd::binary_accumulator::unit_test{

  void accumulating_ones(){

    BinaryAccumulator<Real> acc;
    for(unsigned long j=0; j< (1ul<<6); ++j){
      acc += 1.;
      for(unsigned long u=0; u< size(acc.Stacks); ++u){
	// std::cerr<<std::setprecision(15)<<acc.Stacks[u]<<" ";
      }
      // std::cerr<<"\n";
    }
    acc.MeanValue();
    // for( auto j: acc.SumOfTheFirstTwoTo ){
    //   std::cerr<<"sums_power_of_two = "<<j<<std::endl;
    // }

  }


}//namespace
