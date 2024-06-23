namespace ffd::user_space{

  template<typename G0_t,
	   typename coord_t,
	   typename field_t>
  std::vector<field_t>
  fill_G0_matrix(G0_t const& G0,
		 std::vector<QuantumFieldPair<coord_t>> chis);
  
}//namespace
