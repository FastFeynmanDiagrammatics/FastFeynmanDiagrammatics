

namespace ffd::richardson_extrapolation{

  template<typename Field>
  std::vector<Field>
  RichardsonExtrapolationIntegerPowers(std::vector<Field> sequence){
    std::vector<Field> ret;
    int order = std::size(sequence)-1;
    for(int k=0; k < order; ++k){
      Real pow_2mkm1 = pow(2., -k-1);
      for(int j=0; j < order - k; ++j){
	sequence[j] = (sequence[j+1] - pow_2mkm1*sequence[j])/(1. - pow_2mkm1);
	if(j == order - k - 1){
	  ret.push_back(sequence[j]);
	}
      }
    }
    return ret;
  }

}//namespace
