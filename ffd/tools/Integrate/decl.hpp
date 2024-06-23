namespace ffd::integrate{
  

  enum class Method : int {ClenshawCurtis}; 

  
  template<Method method_name, typename real_function_t>
  Real Integrate(real_function_t FunctionToIntegrate,
		 std::array<Real, 2> IntervalOfIntegration,
		 Real RequestedAbsolutePrecision = 1e-10);
  
	
}//namespace
