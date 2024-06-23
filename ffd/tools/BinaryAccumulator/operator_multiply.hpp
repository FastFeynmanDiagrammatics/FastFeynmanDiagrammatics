

namespace ffd::binary_accumulator{


  template<typename ElementType>
  
  BinaryAccumulator<ElementType>&
  
  operator*=(BinaryAccumulator<ElementType>& accumulator,
	     ElementType element){
    for( auto &x: accumulator.Stacks){
      x *= element;
    }
    for( auto &x: accumulator.SumOfTheFirstTwoTo ){
      x *= element;
    }

    return accumulator;
  }


  
  template<typename ElementType>
  
  BinaryAccumulator<ElementType>&
  
  operator/=(BinaryAccumulator<ElementType>& accumulator,
	     ElementType element){
    return accumulator *= 1./element;
  }

  
  


}//namespace
