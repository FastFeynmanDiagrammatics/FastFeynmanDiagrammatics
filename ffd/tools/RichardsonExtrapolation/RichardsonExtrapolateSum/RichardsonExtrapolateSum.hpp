

namespace ffd::richardson_extrapolation::sums{

  template<typename Field>
  Field
  RichardsonExtrapolateSum(std::function<Field(int)> summand,
			   std::array<const char*, 2> interval,//"-Infinity", "0", "Infinity"
			   [[maybe_unused]] Real RequestedAbsolutePrecision){
    assert( std::string(interval[1]) == "Infinity");
    auto vector_partial_sums = ComputePartialSums(summand, interval[0]);
    
  }
  

}//namespace
