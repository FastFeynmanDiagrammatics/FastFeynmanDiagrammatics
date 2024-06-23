

namespace ffd::binary_accumulator{

  template<typename ElementType>
  
  BinaryAccumulator<ElementType>&
  
  operator+=(BinaryAccumulator<ElementType>& accumulator,
	     ElementType element){
    auto const size_acc = size(accumulator);

    bool SizeIsAPowerOfTwo =
      1ul<<ffd::math_tools::log2_int(size_acc+1)
      == size_acc + 1;
    

    if( SizeIsAPowerOfTwo ){
      accumulator.Stacks.push_back( 0. );
    }

    
    accumulator.Stacks[0] += element;
    for(auto j = 0ul; size_acc & (1ul<<j); ++j){
      accumulator.Stacks[j+1] += accumulator.Stacks[j];
      accumulator.Stacks[j] = 0.;
    }

    
    accumulator.NumElements++;

    if( SizeIsAPowerOfTwo ){
      accumulator.SumOfTheFirstTwoTo.
	push_back(accumulator.
		  Stacks[ffd::math_tools::log2_int(size_acc+1)]);
    }

    
    return accumulator;
  }

  
  
}//namespace
