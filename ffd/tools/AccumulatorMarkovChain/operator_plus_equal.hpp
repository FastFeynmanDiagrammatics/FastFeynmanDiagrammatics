

namespace ffd::accumulator_markov_chain{

  template<typename Field>

  AccumulatorMarkovChain<Field>&
  
  operator+=(AccumulatorMarkovChain<Field>& Acc,
	     Field value){
    unsigned long const max_auto_time = Acc.MaxAutocorrelationTime;


    Acc.Values += value;

    
    for(int j = max_auto_time - 1; j > 0; --j){
      Acc.ValuesHistory[j] = Acc.ValuesHistory[j-1];
    }
    Acc.ValuesHistory[0] = value;

    
    auto value_conj = ffd::complex_conjugate::
      ComplexConjugate( value );
    auto history_min = std::min(max_auto_time,
				size(Acc.Values));
    for( int j: ffd::vector_range::Range(history_min) ){
      Acc.CorrelationFunction[j] += value_conj * Acc.ValuesHistory[j];
    }

    
    return Acc;
  }

  
}//namespace
