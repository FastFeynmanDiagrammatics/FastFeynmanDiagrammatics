

namespace ffd::integral_function{

  template<typename Field>
  class IntegralFunction{
  public:
    std::function<Field(Real)> Integrand;

    Real LowerBound;

    Real RequestedAbsolutePrecision;
    int MinimalOrder = 2;
    Real AddendIncreasingIterativeOrder = 1.;
    Real FactorIncreasingIterativeOrder = 1.51;


    IntegralFunction(std::function<Field(Real)> integrand_,
		     Real lower_bound_,
		     Real requested_precision_ = 1e-10):
      Integrand(integrand_),
      LowerBound(lower_bound_),
      RequestedAbsolutePrecision(requested_precision_){}
      

    Field operator()(Real upper_bound_) const;
    
  };
  

}//namespace
