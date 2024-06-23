

namespace ffd::conformal_mappings{

  template<typename Field>
  class ZinnJustin2{
    std::vector<Field> Coef_w;
    
  public:

    Real CoefficientMapping(int n, int m);
    
    ZinnJustin2(std::vector<Field> const&);

    Complex operator()(Complex) const;
    
    
  };


}//namespace
