namespace ffd::four_point_vertex_hubbard{
	
  template<int order,
	   typename G0_t,
	   typename nilpoly_t =
	   ffd::nilpotent_polynomial::NilpotentPolynomial<Real>>
  
  std::array<nilpoly_t, order*order*order*order>
  connected(G0_t const& G0);


  template<int order,
	   typename G0_t,
	   typename nilpoly_t =
	   ffd::nilpotent_polynomial::NilpotentPolynomial<Real>>
  
  std::array<nilpoly_t, order*order*order*order>
  A_func(G0_t const& G0,
	 std::array<nilpoly_t, order*order*order*order> const& conn);


  
  template<int order,
	   typename G0_t,
	   typename nilpoly_t =
	   ffd::nilpotent_polynomial::NilpotentPolynomial<Real>>
  
  std::array<nilpoly_t, order*order*order*order>
  irreducible(G0_t const& G0,
	      std::array<nilpoly_t, order*order*order*order> const& a_func);

	
}//namespace
