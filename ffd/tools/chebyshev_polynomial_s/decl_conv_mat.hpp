namespace ffd::chebyshev_polynomial_s{

  template<std::size_t n,
	   ffd::phys::statistics stat_zeta = ffd::phys::fermi>
  struct TPoly_conv_mat{
    inline static std::array<Real, n*n*n> conv_mat;
    inline static bool conv_mat_is_computed = false;
  };

}//namespace
