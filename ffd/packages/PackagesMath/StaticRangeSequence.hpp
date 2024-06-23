

namespace ffd::packages_math{

  using IntType = int;
  
  template<IntType...>
  struct StaticRangeSequence;

  template<IntType low, IntType high, IntType... sequence>
  struct StaticRangeSequence<low, high, sequence...>:
    public StaticRangeSequence<low, high-1, high-1, sequence...>{};
  
  template<IntType low, IntType... sequence>
  struct StaticRangeSequence<low, low, sequence...>{};


  template<IntType first, IntType... last>
  std::vector<IntType> PrintSequence(){
    std::vector<IntType> ret;
    ((ret.push_back(first)), ... , (ret.push_back(last)));
    return ret;
  }

  template<IntType low, IntType... sequence>
  std::vector<IntType> PrintSequenceRange(StaticRangeSequence<low, low, sequence...>){
    return PrintSequence<sequence...>();
  }

}//namespace ffd::packages_math
