

namespace ffd::integral_function{

  template<typename Field>
  class IntegralFunction{
  public:
    std::function<Field(Real)> Integrand;

    Real LowerBound;

    Real RequestedAbsolutePrecision = 1e-10;
    int MinimalOrder = 2;
    Real AddendIncreasingIterativeOrder = 1.;
    Real FactorIncreasingIterativeOrder = 1.51;


    IntegralFunction(std::function<Field(Real)> integrand_,
		     Real lower_bound_,
		     Real requested_absolute_precision = 1e-10,
		     int minimal_order = 2,
		     Real addend = 1,
		     Real factor = 1.51):
      Integrand(integrand_), LowerBound(lower_bound_){}
      

    Field operator()(Real upper_bound_) const;
    
  };
  

}//namespace
