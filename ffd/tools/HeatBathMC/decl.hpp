namespace ffd::heat_bath_mc{
  
  using default_tmp_t = ffd::type_traits::default_tmp_t;
  inline constexpr int order_minus_default = 2;
  inline constexpr bool is_translation_invariant_default = false;

  
  template<int max_order,
	   int target_order = max_order-order_minus_default,
	   bool is_translation_invariant = is_translation_invariant_default,
	   int max_calculated_order = max_order,
	   typename nilpotent_polynomial_or_array_t = default_tmp_t,
	   typename seed_function_t = default_tmp_t,
	   typename vector_spacetime_coordinates_t = default_tmp_t,
	   typename origin_spacetime_coordinate_t = default_tmp_t>

  std::pair<std::array<Real, (1<<max_order)>,
    std::array<ffd::type_traits::element_of_container_t<vector_spacetime_coordinates_t>, max_order>>

  HeatBathMC(nilpotent_polynomial_or_array_t const& NilPolynomial_or_array,
	     seed_function_t const& Creates_new_coordinates_and_evaluable,
	     vector_spacetime_coordinates_t const& Vertex_coordinates,
	     origin_spacetime_coordinate_t const& Origin = default_tmp_t());
	
}//namespace
