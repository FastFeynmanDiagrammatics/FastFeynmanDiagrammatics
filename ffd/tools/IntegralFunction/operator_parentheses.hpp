

namespace ffd::integral_function{

  template<typename Field>
  Field IntegralFunction<Field>::operator()(Real upper_bound_) const{
    return
      ffd::integrate_clenshaw_curtis::IntegrateWithClenshawCurtis
      (Integrand,
       {LowerBound, upper_bound_},
       RequestedAbsolutePrecision,
       MinimalOrder,
       AddendIncreasingIterativeOrder,
       FactorIncreasingIterativeOrder);
  }

  
}//namespace
