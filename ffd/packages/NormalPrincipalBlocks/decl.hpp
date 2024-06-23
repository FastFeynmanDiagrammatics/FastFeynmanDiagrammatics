namespace ffd::normal_principal_blocks{
  
  using namespace ffd::user_space;
  
  
  template<typename coord_t, typename field_t>
  std::vector<std::pair<std::vector<QuantumFieldPair<coord_t>>,
			std::array<std::vector<int>, 2>>>
  Blocks(QuantumFieldSum<coord_t, field_t> const& interaction,
	 QuantumFieldSum<coord_t, field_t> const& external_vertices,
	 QuantumFieldSum<coord_t, field_t> const& action);
  
  
  
  template<typename timer_t,
    typename nilpoly_t =
	   ffd::nilpotent_polynomial::NilpotentPolynomial<Real>,
	   typename coord_t=int, typename field_t=int,
	   typename G0_t=int,
	   typename alpha_shift_t=int
	   >
  std::array<nilpoly_t, 2>
  PM(timer_t& timer_v,
     std::vector<QuantumFieldSum<coord_t, field_t>> const& interaction,
     QuantumFieldSum<coord_t, field_t> const& external_vertices,
     QuantumFieldSum<coord_t, field_t> const& action,
     G0_t const& G0,
     alpha_shift_t const& alpha_shift = 0);
  
  
  
}//namespace
