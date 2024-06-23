namespace ffd::accumulator_markov_chain{

  template<typename Field = Real>
  
  struct AccumulatorMarkovChain{

    using bin_acc_t = ffd::binary_accumulator::
      BinaryAccumulator<Field>;

    
    bin_acc_t Values;

    
    int MaxAutocorrelationTime,
      TailAutocorrelationTime = 8;
    std::vector<Field> ValuesHistory;
    std::vector<bin_acc_t> CorrelationFunction;
    
    
    
    AccumulatorMarkovChain(Real MaxAutocorrelationTime_ = 100):
      MaxAutocorrelationTime(MaxAutocorrelationTime_),
      ValuesHistory(MaxAutocorrelationTime_, 0.),
      CorrelationFunction(MaxAutocorrelationTime_) {}


    Real
    Variance() const;


    Real
    AutocorrelationFunction(int time) const;
    
    
    Real
    AutocorrelationTime() const;
    
    
    
    std::pair<Field, Real>
    MeanValueError() const;
    
  };


  template<typename Field>
  auto size(AccumulatorMarkovChain<Field> const& acc){
    return size(acc.Values);
  }

  
}//namespace
