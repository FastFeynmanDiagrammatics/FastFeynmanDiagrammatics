

namespace ffd::richardson_extrapolation{

  std::vector<Complex>
  RichardsonExtrapolationIntegerPowersPhase(std::vector<Complex> sequence, Real Phase){
    std::vector<Complex> ret;
    int order = std::size(sequence)-1;
    for(int k=0; k < order; ++k){
      for(int j=0; j < order - k; ++j){
	Complex pow_2mkm1_phase = pow(2., -k-1)*exp(Complex(0, Phase*(1<<(j+k))));
	sequence[j] = (sequence[j+1] - pow_2mkm1_phase*sequence[j])/(1. - pow_2mkm1_phase);
	if(j == order - k - 1){
	  ret.push_back(sequence[j]);
	}
      }
    }
    return ret;
  }

}//namespace
