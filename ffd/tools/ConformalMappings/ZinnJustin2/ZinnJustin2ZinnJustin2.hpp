

namespace ffd::conformal_mappings{

  template<typename Field>
  ZinnJustin2<Field>::ZinnJustin2(std::vector<Field> const& coef_z){
    Coef_w = std::vector<Field>(std::size(coef_z));
    Coef_w[0] = coef_z[0];
    for(std::size_t m=1; m < std::size(Coef_w); ++m){
      Coef_w[m] = 0;
      for(std::size_t n=m; n >= 1; --n){
	Coef_w[m] += coef_z[n]*CoefficientMapping(n, m);
      }
    }
  }
  

}//namespace
