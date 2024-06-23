namespace ffd::fft_matsubara_3d {

template <std::floating_point float_t = double>
auto OmegaToTau_FFF(ffd::vector3d<std::complex<float_t>> M, float_t Beta) {
  using complex_t = std::complex<float_t>;
  auto const n_omega = size_x(M);
  assert((size_x(M) == size_y(M)));
  assert((size_x(M) == size_z(M)));
  ffd::vector3d<complex_t> M_inv(n_omega, n_omega, n_omega);
  for (long j = 0; j < size(M); ++j) {
    M_inv[j] = M[n_omega * n_omega * n_omega - j - 1];
    // M(u, k, l) = M(n-u, n-k, n-l)
    // M[j] = M[n^3 - j - 1]
  }  //  for j in range(0, F.size())
  auto T = M;
  const char* err_str = NULL;
  bool res;
  res = simple_fft::FFT(M_inv, T, n_omega, n_omega, n_omega, err_str);
  auto n_omega_eff = n_omega / 2;
  auto tau_k = [Beta, n_omega](auto k) { return k * Beta / n_omega; };
  for (ulong z = 0; z < n_omega; ++z) {
    for (ulong y = 0; y < n_omega; ++y) {
      for (ulong x = 0; x < n_omega; ++x) {
        T(x, y, z) *=
            std::exp(-2 * M_PI * (n_omega_eff + .5) *
                     (tau_k(x) + tau_k(y) + tau_k(z)) * complex_t(0., 1.)) /
            std::pow(Beta, 3);
      }  // for x in range(0, n_omega)
    }    // for y in range(0, n_omega)
  }      // for z in range(0, n_omega)
  return T;
}

template <std::floating_point float_t = double>
auto OmegaToTau_FFB(ffd::vector3d<std::complex<float_t>> M, float_t Beta) {
  using complex_t = std::complex<float_t>;
  auto const N_F_tot = size_x(M);
  assert((size_x(M) == size_y(M)));
  auto const N_B_tot = size_z(M);
  ffd::vector3d<complex_t> M_inv(N_F_tot, N_F_tot, N_B_tot);
  for (ulong j = 0; j < N_F_tot; ++j) {
    for (ulong k = 0; k < N_F_tot; ++k) {
      for (ulong l = 0; l < N_B_tot; ++l) {
        M_inv(j, k, l) = M(N_F_tot - 1 - j, N_F_tot - 1 - k, N_B_tot - 1 - l);
      }  // for l in range(0, N_B_tot)
    }    // for k in range(0, N_F_tot)
  }      // for j in range(0, N_F_tot)
  auto T = M;
  const char* err_str = NULL;
  bool res;
  res = simple_fft::FFT(M_inv, T, N_F_tot, N_F_tot, N_B_tot, err_str);
  auto const N_F = N_F_tot / 2;
  auto const N_B = N_B_tot / 2 - 1;
  auto tau_F = [Beta, N_F_tot](auto k) { return k * Beta / N_F_tot; };
  auto tau_B = [Beta, N_B_tot](auto k) { return k * Beta / N_B_tot; };
  for (ulong z = 0; z < N_B_tot; ++z) {
    for (ulong y = 0; y < N_F_tot; ++y) {
      for (ulong x = 0; x < N_F_tot; ++x) {
        T(x, y, z) *=
            std::exp(-2 * M_PI *
                     ((N_F + .5) * (tau_F(x) + tau_F(y)) + N_B * tau_B(z)) *
                     complex_t(0., 1.)) /
            std::pow(Beta, 3);
      }  // for x in range(0, n_omega)
    }    // for y in range(0, n_omega)
  }      // for z in range(0, n_omega)
  return T;
}

}  // namespace ffd::fft_matsubara_3d
