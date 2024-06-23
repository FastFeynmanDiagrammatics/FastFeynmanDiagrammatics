namespace ffd::accumulator_markov_chain{

  template<typename Field>

  Real

  AccumulatorMarkovChain<Field>::
  AutocorrelationTime() const{
    Real tau_int = 0;

    
    auto const num_samples = size( Values );
    

    for(int j=1; j < MaxAutocorrelationTime-TailAutocorrelationTime; ++j){
      tau_int += AutocorrelationFunction(j)*(1-j*1./num_samples);
    }

    Real tau_tail = 0;
    for(int j = MaxAutocorrelationTime-TailAutocorrelationTime;
	j < MaxAutocorrelationTime; ++j){
      tau_tail += AutocorrelationFunction(j)*(1-j*1./num_samples);
    }
    Real alpha_approx = std::pow( std::abs(AutocorrelationFunction(MaxAutocorrelationTime)),
			     1./(MaxAutocorrelationTime-1) );
    
    return tau_int + tau_tail/(1-std::pow(alpha_approx,
					  TailAutocorrelationTime));
  }

}//namespace
