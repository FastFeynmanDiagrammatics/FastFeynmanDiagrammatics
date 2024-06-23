namespace ffd::fft_matsubara_3d::unit_test {

void simple() {
  using cpx_t = std::complex<double>;
  ffd::vector3d<cpx_t> v(4, 4, 8, 1.);
  auto v_fourier = OmegaToTau_FFB(v, 1.5);
  for (ulong x = 0; x < size_x(v_fourier); ++x) {
    for (ulong y = 0; y < size_y(v_fourier); ++y) {
      for (ulong z = 0; z < size_z(v_fourier); ++z) {
        std::cerr << x << ' ' << y << ' ' << z << ' ' << v_fourier(x, y, z)
                  << '\n';
      }  // for z in range(0, size_z(v_fourier))
    }    // for y in range(0, size_y(v_fourier))
  }      // for x in range(0, size_x(v_fourier))
}

}  // namespace ffd::fft_matsubara_3d::unit_test
