namespace ffd::binary_accumulator{

  template<typename element_t = Real>
  struct BinaryAccumulator{

    
    unsigned long NumElements;
    std::vector<element_t> Stacks;
    std::vector<element_t> SumOfTheFirstTwoTo;

    BinaryAccumulator(): NumElements{0},
			 Stacks{},
			 SumOfTheFirstTwoTo{} {}

    
    element_t
    MeanValue() const;
    
  };

    
  template<typename T>
  unsigned long size(BinaryAccumulator<T> const& X){
    return X.NumElements;
  }


  BinaryAccumulator() -> BinaryAccumulator<Real>;
  
	
}//namespace
