namespace ffd::pade {

template <class field = Real, class vector0_t = int>
auto eval(int p_max, int q_max, vector0_t const& coef) {
  assert((size(coef) - 1 > p_max + q_max));
  auto const lin_size = p_max + 1 + q_max;
  ffd::vector2d<field> M(lin_size, lin_size, field(0.));
  for (ulong j = 0; j < p_max + 1; ++j) {
    M(j, j) = 1.;
  }  // for j in range(0, p_max+1)
  for (ulong j = 0; j < lin_size; ++j) {
    for (ulong l = 1; l <= std::min<int>(j, q_max); ++l) {
      M(p_max + l, j) = -coef[j - l];
    }  // for l in range(0, j)
  }

  // for (ulong j = 0; j < lin_size; ++j) {
  //   for (ulong k = 0; k < lin_size; ++k) {
  //     std::cerr << M(j, k) << " ";
  //   }  // for k in range(0, lin_size)
  //   std::cerr << "\n";
  // }  // for j in range(0, lin_size)

  auto [I, det] = ffd::inverse_matrix_gauss::Inverse_and_Determinant(M);

  // for (ulong j = 0; j < lin_size; ++j) {
  //   for (ulong k = 0; k < lin_size; ++k) {
  //     std::cerr << I[j * lin_size + k] << " ";
  //   }  // for k in range(0, lin_size)
  //   std::cerr << "\n";
  // }  // for j in range(0, lin_size)

  ffd::vector1d<field> p(p_max + 1), q(q_max);
  for (ulong j = 0; j <= p_max; ++j) {
    p[j] = 0.;
    for (ulong k = 0; k < lin_size; ++k) {
      p[j] += I[j * lin_size + k] * coef[k];
    }  // for k in range(0, lin_size)
  }    // for j in range(0, lin_size)
  for (ulong j = 0; j < q_max; ++j) {
    q[j] = 0.;
    for (ulong k = 0; k < lin_size; ++k) {
      q[j] += I[(p_max + 1 + j) * lin_size + k] * coef[k];
    }  // for k in range(0, lin_size)
  }    // for j in range(0, lin_size)

  return std::make_pair(p, q);
}

template <class field, class vector0_t, class vector1_t>
field eval_rational_function(vector0_t const& p, vector1_t const& q, field z) {
  field P = 0., Q = 1., pow_z = 1.;
  for (ulong j = 0; j < size(p); ++j, pow_z *= z) {
    P += p[j] * pow_z;
  }  // for j in range(0, size(p))
  pow_z = z;
  for (ulong j = 0; j < size(q); ++j, pow_z *= z) {
    Q += q[j] * pow_z;
  }  // for j in range(0, size(p))
  return P / Q;
}

}  // namespace ffd::pade
