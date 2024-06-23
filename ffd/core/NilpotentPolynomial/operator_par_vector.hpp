namespace ffd::nilpotent_polynomial{

  template <typename Field, SetT SetType>
  
  NilpotentPolynomial<Field, SetType>
  NilpotentPolynomial<Field, SetType>::
  operator()(std::vector<Field> const& int_vec) const{
    auto ret = *this;
    for(BinaryInt S=0; S<ret.size(); ++S){
      auto v = ffd::set_theory::VectorOfBinaryDigitsOf(S);
      Field factor = 1.;
      for(auto const& x: v) factor *= int_vec[x];
      ret[S] *= factor;
    }
    return ret;
  }

}//namespace
