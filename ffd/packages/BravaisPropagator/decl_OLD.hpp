namespace ffd::bravais_propagator{

  class BravaisPropagator{
  public:

    using chebyshev_t = ffd::chebyshev_polynomial::ChebyshevPolynomial<Real>;
    
  private:
    Real const precision_G0_default = 1e-8;
    std::vector<chebyshev_t> G0;
    Real beta_;
    std::vector<int> L;
      
  public:


    template<typename action_t,
	     typename time_t,
	     typename lattice_t>
    BravaisPropagator(action_t const& S0,
		      time_t,
		      lattice_t,
		      Real precision_G0 = -1);


    template<typename coord_t>
    [[nodiscard]] Real
    operator()(coord_t const& X0, coord_t const& X1, int spin = 0) const;
    
  };
  
}//namespace
