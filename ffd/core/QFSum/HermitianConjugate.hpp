namespace ffd::user_space{

  template<class coord_d,
	   class field_d>
  auto
  HermitianConjugate(QFSum<coord_d, field_d> const& x){
    auto ret = Bar(x);
    if constexpr(std::is_same_v<field_d, Complex>){
    for(uint j=0; j<size(ret.fields); ++j){
      ret.coef[j] = std::conj(ret.coef[j]);
    }
    }
    return ret;
  }

}//namespace
