namespace ffd::find_root{

  template<typename function_t>
  
  std::optional<Real>
  
  FindRealRoot(function_t real_function,
	       Real starting_point,
	       Real AbsolutePrecision = parameters::absolute_precision,
	       std::array<std::optional<Real>, 2> domain =

	       std::array{std::optional<Real>{}, std::optional<Real>{}});


  

  template<typename function_t>
  
  std::optional<Real>
  
  FindRealRoot(function_t real_function,
	       std::array<Real, 2> bracket,
	       Real AbsolutePrecision = parameters::absolute_precision);

  
  template<int d,
	   typename function_t>

  std::optional<std::array<Real, d>>

  FindRealRoots(function_t real_vector_function,
		std::array<Real, d> starting_points,
		Real AbsolutePrecision = parameters::absolute_precision);

  
	
}//namespace
