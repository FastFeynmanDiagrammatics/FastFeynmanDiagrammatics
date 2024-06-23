

namespace ffd::binary_accumulator{

  template<typename ElementType>

  ElementType
  
  BinaryAccumulator<ElementType>::
  MeanValue() const{
    ElementType ret = 0.;

    
    for(auto const& stack: Stacks){
      ret += stack/Real(NumElements);
    }


    return ret;
  }

}//namespace
