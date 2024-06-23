

namespace ffd::richardson_extrapolation::sums{

  template<typename Field>
  std::vector<Field> ComputePartialSums(std::function<Field(int)> summand,
					int log2_maximal_int,
					const char* lower_bound = "-Infinity"){
    const bool IsLowerBoundInfinity = std::string(lower_bound) == "-Infinity";
    const bool IsLowerBoundOne = std::string(lower_bound) == "1";
    assert(IsLowerBoundOne || IsLowerBoundInfinity || std::string(lower_bound) == "0");
    std::vector<Field> ret(1<<log2_maximal_int);
    if(IsLowerBoundOne){
      ret[0] = 0;
    }else{
      ret[0] = summand(0);
    }
    for(int u = 1; u < log2_maximal_int; ++u){
      ret[u] = ret[u-1];
      for(long j = 1<<(u-1); j < 1<<u; ++j){
	ret[u] += summand(j);
	if(IsLowerBoundInfinity){
	  ret[u] += summand(-j);
	}
      }
    }
    return ret;
  }

}//namespace
