namespace ffd::conf_map {

template <class cmap_f>
auto root(cmap_f const& cmap, Real U) {
  Real wm = 0, wM = 1, wmean;
  do {
    wmean = .5 * (wm + wM);
    if (auto U_now = cmap(wmean); U_now < U) {
      wm = wmean;
    } else {
      wM = wmean;
    }
  } while (wM - wm > 1e-12);
  return wmean;
}  // namespace ffd::conf_map

}  // namespace ffd::conf_map
