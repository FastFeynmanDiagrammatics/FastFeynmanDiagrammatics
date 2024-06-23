namespace ffd::canonical_sigma {

template <std::size_t n, class field_d = Real>

void normal_pow_alpha_r(array2d<field_d, n, n> const& __restrict__ G0,
                        vector2d<field_d>& __restrict__ pow_alpha) {
  auto const PM = ffd::principal_minors::PrincipalMinorsBFS(G0);
  pow_alpha.fill(0);

  return;
}

}  // namespace ffd::canonical_sigma
