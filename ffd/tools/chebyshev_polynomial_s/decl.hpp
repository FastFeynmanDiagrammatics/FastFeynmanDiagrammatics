namespace ffd::chebyshev_polynomial_s{

  template<std::size_t n,
	   typename field_t = Real,
	   ffd::phys::statistics stat_zeta = ffd::phys::fermi,
	   template<typename, std::size_t> typename array_t = std::array>
  struct TPoly{
    array_t<field_t, n> coef;
    Real bounds_average, one_over_half_diff;
    inline static bool cosinus_are_computed = false;
    inline static std::array<Real, n*n> cosinus_transform;
    
    TPoly(): coef() {}

    TPoly(std::array<Real, 2> lower_upper_limit){coef.fill(0.);
      bounds_average = .5*(lower_upper_limit[0]+lower_upper_limit[1]);
      one_over_half_diff = 2./(-lower_upper_limit[0]+lower_upper_limit[1]);}

    TPoly(std::array<field_t, n> const& values_at_chebyshev_nodes,
	  std::array<Real, 2> lower_upper_limit);

    
    template<typename function_t>
    TPoly(function_t const& f,
	  std::array<Real, 2> lower_upper_limit);


    [[nodiscard]] field_t operator()(Real const x) const noexcept;

    template<class stream_t>
    void FromStream(stream_t&);
    
  };


  template<std::size_t n, typename field_t>
  constexpr std::size_t
  size(TPoly<n, field_t> const&){ return n;}
	
}//namespace
