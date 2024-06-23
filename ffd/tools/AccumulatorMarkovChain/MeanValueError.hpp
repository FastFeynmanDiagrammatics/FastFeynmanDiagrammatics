namespace ffd::accumulator_markov_chain{

  
  template<typename Field>

  std::pair<Field, Real>

  AccumulatorMarkovChain<Field>::
  MeanValueError() const{
    Real mean_value = Values.MeanValue();


    Real effective_variance = Variance()*(1 + 2*AutocorrelationTime() );
    unsigned long samples = size( Values );
    
    
    Real error = 
      std::sqrt(std::abs(effective_variance/samples));


    return std::pair{mean_value, error};    
  }


}//namespace
