namespace ffd::conf_map {

template <class field, auto n, template <class, auto> class P, class cmap_f>
auto map_poly(cmap_f const& cmap, P<field, n> const& p) {
  auto w = P<field, n>(Real(0));
  w[1] = 1;
  auto const cw = cmap(w);
  auto pow_cw = cw;
  auto ret = P<field, n>{p[0]};
  for (ulong j = 1; j < p.order + 1; ++j) {
    ret += p[j] * pow_cw;
    pow_cw *= cw;
  }  // for j in range(0, p.order+1)
  return ret;
}

}  // namespace ffd::conf_map
