namespace ffd::accumulator_markov_chain{

  template<typename Field>

  Real

  AccumulatorMarkovChain<Field>::
  Variance() const{
    Real variance = ffd::real_complex_tools::
      real(CorrelationFunction[0].MeanValue());
    
    
    variance -= std::pow( std::abs(Values.MeanValue()), 2);


    return variance >= 0. ? variance : 0.;
  }

}//namespace
