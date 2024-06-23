namespace ffd::heat_bath_mc{

  template<uint n, uint k>
  auto
  compute_bin_coef(){
    std::array<Real, k+1> ret;
    for(uint j=0; j<k+1; ++j){
      ret[j] = 1./ffd::core_math::BinomialCoefficient(n, j);
    }
    return ret;
  }

}//namespace
