namespace ffd::normal_propagator_g {

template <int d,
          int n_atoms = 1,
          typename field_d = Real,
          typename Eig_k_blocks_d = int,
          typename k_array_d = int,
          typename block_oracle_f = int>
auto exact_propagator_tau_0T(
    std::array<int, d> const& L,
    Eig_k_blocks_d const& Eig_k_blocks,
    k_array_d const& k_array,
    block_oracle_f const& b_oracle,
    uint u,  // block number
    std::size_t p0,
    std::size_t a0,  // p0: particle in block u, a0: atom
    std::size_t p1,
    std::size_t a1,  // p1: particle in block u, a1: atom
    std::array<int, d>
        x0) {  // relative brillouin coordinate of the zeroth particle

  using namespace ffd::user_space;
  using ffd::core_math::Pi;

  // building H_{k,\sigma}
  auto blocks = b_oracle();
  auto ret_index = ffd::index_lattice::g<d, n_atoms>(L, b_oracle);
  auto index_f = ret_index.first;

  std::size_t const size_block = blocks.split[u + 1] - blocks.split[u];
  std::size_t const n_comp = size_block * n_atoms;

  auto const s0 = p0 + size_block * a0;
  auto const s1 = p1 + size_block * a1;
  auto const& Eig_k = Eig_k_blocks[u];
  Real one_L_prod = 1;
  for (int j = 0; j < d; ++j)
    one_L_prod /= L[j];

  auto G0_tau = [s0, s1, x0, L, one_L_prod, k_array,
                 Eig_k](Real tau) -> field_d {
    Real const tau_sign = 1 - 2 * (tau <= 0);
    //	std::cerr << "tau_sign = " << tau_sign << '\n';
    Complex ret = 0.;
    for (std::size_t k_ind = 0; k_ind < size(Eig_k); ++k_ind) {
      auto k = k_array[k_ind];
      auto e_k = Eig_k[k_ind];
      auto E_k = e_k.first;
      auto v_k = e_k.second;
      Complex ret_k = 0.;
      for (uint j = 0; j < size(E_k); ++j) {
        //	    std::cerr << "E_k_j = " << j << ' ' << E_k[j] << ' ' <<
        //tau*E_k[j] << '\n';
        Real const tau_E = tau * E_k[j];
        //////ONLY WORKS FOR FERMIONS right now
        if (tau_E >= 0) {
          ret_k += tau_sign * v_k[j][s0] * conj(v_k[j][s1]) * exp(-tau_E);
        }
        // else if( std::abs(E_k[j]) == 0 ){
        //   ret_k += tau_sign*v_k[j][s0]*conj(v_k[j][s1]);
        // }
      }

      Real k_x = 0;
      for (int j = 0; j < d; ++j) {
        k_x += (k[j] * x0[j] * 2.) * Pi / L[j];
      }

      ret += std::exp(Complex(0., k_x)) * ret_k;
    }

    if constexpr (std::is_same_v<field_d, Real>) {
      return std::real(ret * one_L_prod);
    } else {
      return ret * one_L_prod;
    }
  };

  return G0_tau;
}

}  // namespace ffd::normal_propagator_g
