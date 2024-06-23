namespace ffd::accumulator_markov_chain{

  template<typename Field>

  Real

  AccumulatorMarkovChain<Field>::
  AutocorrelationFunction(int time) const{

    
    Real const variance = Variance();
    Real const mean_square = std::pow( std::abs(Values.MeanValue()), 2);


    Real corr_func = size(CorrelationFunction) > unsigned(time) ? ffd::real_complex_tools::
      real( CorrelationFunction[time].MeanValue() ) : 0.;


    corr_func -= mean_square;
    corr_func /= (variance > 0 ? variance : 1e-200);

    
    return corr_func;
  }

}//namespace
