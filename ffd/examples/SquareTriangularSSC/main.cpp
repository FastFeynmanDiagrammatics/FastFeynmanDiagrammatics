#include <ffd.hpp>
#include "constants.hpp"
#include "parameters.hpp"

using namespace ffd::user_space;
namespace fs = std::filesystem;

int main() {
  auto const [process_ID, thread_ID] = [] {
    std::stringstream ss_tid, ss_pid;
    ss_pid << pid_counter;
    ss_tid << ss_pid.str() << "_" << RandomString(10) << ".ds";
    {
      std::ofstream pid_file("parameters/pid_counter");
      pid_file << pid_counter + 1;
    }
    fs::remove("kill");
    fs::create_directory("out");
    fs::current_path("out");
    fs::remove("kill");
    fs::create_directory(ss_pid.str());
    fs::current_path(ss_pid.str());
    fs::remove("kill");
    fs::create_directory("parameters");
    fs::copy("../../parameters", "parameters",
             fs::copy_options::overwrite_existing);
    fs::copy("../../parameters.hpp", ".", fs::copy_options::overwrite_existing);
    fs::copy("../../main.cpp", ".", fs::copy_options::overwrite_existing);
    fs::create_directory("lambda_norm");
    fs::current_path("..");
    return std::make_pair(ss_pid.str(), ss_tid.str());
  }();

  SafeArray<Timer, (1ul << 8)> timers;
  constexpr int d = 2;
  int normalization_phase = 1;
  using accum_t = ffd::blocking_method::block_t<1024, 16, Real>;
  using poly_t = ffd::truncated_polynomial::P<actual_order, Real>;

  auto const [H0, mu0] = [] {
    auto x = LatticeCoordinates_g<d>();
    auto hterm = [x](int x0, int x1) {
      return Bar(psi('u')(x(x0, x1))) * psi('u')(x(0, 0));
    };
    auto H_hop_spin_t = [hterm] {
      auto ret = hterm(1, 0) + hterm(0, 1);
      if (geometry == triangle)
        ret += hterm(-1, 1);
      return ret;
    }();
    auto H_hop_t = -t_hopping * H_hop_spin_t;
    H_hop_t += HermitianConjugate(H_hop_t);

    auto H_hop_spin_tp = [hterm] {
      auto ret = hterm(1, 1);
      if (geometry == square)
        ret += hterm(-1, 1);
      else if (geometry == triangle) {
        ret += hterm(2, 0) + hterm(0, 2) + hterm(-1, 2) + hterm(-2, 2) +
               hterm(2, 1);
      }
      return ret;
    }();
    auto H_hop_tp = -tprime_hopping * H_hop_spin_tp;
    H_hop_tp += HermitianConjugate(H_hop_tp);

    auto H_t = H_hop_t + H_hop_tp;

    auto H0_mu = hterm(0, 0);
    Real mu_min = -10, mu_max = 10, mu = 0, n_0_mu = 0;
    while (mu_max - mu_min > precision_mu) {
      mu = .5 * (mu_min + mu_max);
      auto H_0 = H_t - mu * H0_mu;
      auto [Eig_k_blocks, k_array, block_oracle] =
          ffd::normal_propagator_g::matrix_diag<d>({Lx, Ly}, H_0);
      auto G0_tau_up = ffd::normal_propagator_g::exact_propagator_tau<d>(
          Beta, {Lx, Ly}, Eig_k_blocks, k_array, block_oracle,
          0,        // correlation block number
          0,        // particle index of first field
          0,        // atomic position of first field
          0,        // particle index of second field
          0,        // atomic position of second field
          {0, 0});  // brillouin vector (difference in position)
      n_0_mu = 2 * G0_tau_up(Beta);
      n_0_mu < n0 ? mu_min = mu : mu_max = mu;
    }  // while
    if (verbose_stdout >= 1) {
      std::cout << std::setprecision(10) << "n0(mu=" << mu << ") = " << n_0_mu
                << "\n";
    }
    return std::make_pair(H_t - mu * H0_mu, mu);
  }();

  auto const G0 = [H0 = H0]() {
    auto temp_G0 =
        ffd::normal_propagator_g::g<d>(Beta, {Lx, Ly}, H0, precision_G0);
    return temp_G0.first;
  }();

  auto const [proposer, k_bravais] = [] {
    auto const bravais = [] {
      std::array<std::array<Real, d>, d> ret;
      if (geometry == square) {
        ret[0] = {1., 0.};
        ret[1] = {0., 1.};
      } else if (geometry == triangle) {
        ret[0] = {1., 0.};
        ret[1] = {.5, .5 * std::sqrt(3.)};
      }
      return ret;
    }();
    auto const k_bravais = [b = bravais] {
      std::array<std::array<Real, d>, d> ret;
      auto det = b[0][0] * b[1][1] - b[0][1] * b[1][0];
      ret[0] = {b[1][1] / det, -b[1][0] / det};
      ret[1] = {-b[0][1] / det, b[0][0] / det};
      return ret;
    }();
    auto const Realspace = ffd::lattice_g::Realspace_g<d>({Lx, Ly}, bravais);
    auto const Distance = ffd::lattice_g::Distance_g<d>(Realspace);
    return std::make_pair(
        ffd::itime_lattice_proposer::g(
            ffd::itime_proposer::log_g(Beta, time_scale),
            ffd::lattice_proposer::g<d>(
                {Lx, Ly},
                [](auto x) {
                  auto y = x / space_scale;
                  y *= y;
                  return 1. / (1. + .5 * y + .125 * y * y  //+ y * y * y / 48.
                              );
                },
                Distance)),
        k_bravais);
  }();

  auto X = Spacetime_g<d>();
  auto X_i = [X] {
    std::array<std::decay_t<decltype(X(0., 0, 0))>, order> X_i;
    for (uint j = 0; j < order; ++j)
      X_i[j] = X(Beta * Proba(), 0, 0);

    X_i[0] = X(Beta * 0.1, 0, 0);
    X_i[1] = X(Beta * 0.9, 0, 0);
    X_i[2] = X(Beta * 0.3, 0, 0);
    X_i[3] = X(Beta * 0.7, 0, 0);
    X_i[4] = X(Beta * 0.5, 0, 0);
    X_i[5] = X(Beta * 0.5, 0, 0);
    return X_i;
  }();

  // X_i[0] = X(Beta * 0.5, 0, 0);
  // auto psi_u = psi_f('u');
  // auto bpsi_u = Bar(psi_f('u'));
  // for (uint x = 0; x <= 2 * Lx; ++x) {
  //   X_i[1] = X(Beta * 0.6, x, 0);
  //   X_i[2] = X(Beta * 0.6, 0, x);
  //   std::cout << x << " "
  //             << G0(std::make_pair(psi_u, X_i[0]),
  //                   std::make_pair(bpsi_u, X_i[1]))
  //             << " "
  //             << G0(std::make_pair(psi_u, X_i[0]),
  //                   std::make_pair(bpsi_u, X_i[2]))
  //             << std::endl;
  // }
  // exit(1);

  auto lattice_equivalent = [](std::array<int, 2> n0,
                               std::array<int, 2> n1) -> bool {
    const std::array<int, 2> L{Lx, Ly};
    auto are_equivalent = [L](auto n2, auto n3) -> bool {
      bool ret = true;
      for (ulong s = 0; s < 2; ++s) {
        int diff = n2[s] - n3[s];
        while (diff >= L[s])
          diff -= L[s];
        while (diff < 0)
          diff += L[s];
        ret = ret && diff == 0;
      }  // for s in range(0, 2)
      return ret;
    };
    auto swap_k = [](auto& n) {
      auto t = n[0];
      n[0] = n[1];
      n[1] = t;
    };
    auto pi_rotation = [](auto& n) {
      n[0] = -n[0];
      n[1] = -n[1];
    };
    auto mirror_x = [](auto& n) { n[0] = -n[0]; };
    auto mirror_y = [](auto& n) { n[1] = -n[1]; };
    auto pi_2_rotation = [](auto& n) {
      auto t = n[0];
      n[0] = n[1];
      n[1] = -t;
    };
    auto pi_3_rotation = [](auto& n) {
      auto t = n[0];
      n[0] = n[1];
      n[1] -= t;
    };
    if (Lx == Ly) {
      for (ulong s = 0; s < 2; ++s) {
        swap_k(n0);
        if (geometry == square) {
          for (ulong j = 0; j < 4; ++j) {
            pi_2_rotation(n0);
            if (are_equivalent(n0, n1))
              return true;
          }  // for j in range(0, 4)
        } else if (geometry == triangle) {
          for (ulong j = 0; j < 6; ++j) {
            pi_3_rotation(n0);
            if (are_equivalent(n0, n1))
              return true;
          }  // for j in range(0, 4)
        }
      }  // for s in range(0, 2)
    } else {
      if (geometry == square) {
        mirror_x(n0);
        if (are_equivalent(n0, n1))
          return true;
        mirror_y(n0);
        if (are_equivalent(n0, n1))
          return true;
        mirror_x(n0);
        if (are_equivalent(n0, n1))
          return true;
        return false;
      } else {
        std::cerr << "Wrong geometry settings " << std::endl;
        exit(1);
      }
    }

    return false;
  };

  auto const k_vec = [lattice_equivalent] {
    std::vector<std::array<int, 2>> vector_k;
    for (int fac = 1; fac <= Lx; ++fac) {
      auto fac_x = fac;
      auto fac_y = fac;
      if (Lx < Ly) {
        fac_x = 1;
      } else if (Lx > Ly) {
        fac_y = 1;
      }
      if ((Lx % fac_x == 0) && (Ly % fac_y == 0)) {
        vector_k.resize(0);
        for (int ny = 0; ny < Ly; ny += fac_y) {
          for (int nx = 0; nx < Lx; nx += fac_x) {
            std::array<int, 2> candidate{nx, ny};
            bool not_inside = true;
            for (auto v : vector_k)
              not_inside = not_inside && !lattice_equivalent(candidate, v);
            if (not_inside) {
              vector_k.push_back(candidate);
            }
          }
        }
        if (size(vector_k) <= n_k_points_max)
          break;
      }
    }
    std::vector<std::array<Real, 2>> ret;
    for (ulong j = 0; j < size(vector_k); ++j) {
      ret.push_back(
          {2 * vector_k[j][0] * M_PI / Lx, 2 * vector_k[j][1] * M_PI / Ly});
    }
    if (verbose_stdout >= 1)
      std::cout << "size_k_vec = " << size(ret) << "\n";
    return ret;
  }();

  // for (uint j = 0; j < size(k_vec); ++j) {
  //   std::cout << j << " " << k_vec[j][0] << " " << k_vec[j][1] << std::endl;
  // }
  // exit(1);

  auto const size_k = size(k_vec);
  auto const size_w = omega_end - omega_begin;
  auto const size_kw = size_k * size_w;

  ffd::vector4d<Complex> exp_K(order, order, size_w, size_k, Real(0));
  auto exp_K_r = [&normalization_phase, &k_vec, &exp_K, &X_i]() mutable {
    if (  // (t_MC%accumulate_every) != 0 ||
        normalization_phase != 0)
      return;
    ffd::array2d<Complex, order, omega_end - omega_begin> exp_om;
    for (ulong j = 0; j < size(exp_om); ++j) {
      auto const [u, om] = exp_om.indexes(j);
      exp_om[j] = std::exp(
          Complex(0., -X_i[u].first * (om + omega_begin) * 2. * M_PI / Beta));
    }

    auto exp_k = ffd::vector2d<Complex>{order, size(k_vec)};
    for (ulong j = 0; j < size(exp_k); ++j) {
      auto const [u, k] = exp_k.indexes(j);
      auto const [kx, ky] = k_vec[k];
      auto const [x, y] = X_i[u].second;
      exp_k[j] = std::exp(Complex(0., x * kx + y * ky));
    }

    auto exp_om_k =
        ffd::vector3d<Complex>{order, omega_end - omega_begin, size(k_vec)};
    for (ulong j = 0; j < size(exp_om_k); ++j) {
      auto const [u, om, k] = exp_om_k.indexes(j);
      exp_om_k[j] = exp_k(u, k) * exp_om(u, om);
    }

    for (ulong j = 0; j < size(exp_K); ++j) {
      auto const [u, v, om, k] = exp_K.indexes(j);
      if (u != v) {
        exp_K[j] = exp_om_k(u, om, k) * conj(exp_om_k(v, om, k));
      }
    }
    return;
  };

  auto constexpr p_max = [] {
    std::array<ulong, actual_order + 2> p_max = {0};
    p_max.fill(0);
    for (int j = 1; j <= actual_order; ++j) {
      p_max[j] = std::min<int>(j, actual_order - j);
    }
    // p_max.fill(0);
    return p_max;
  }();

  auto constexpr one_factorial = [] {
    std::array<Real, actual_order + 1> one_factorial;
    for (ulong u = 0; u <= actual_order; ++u)
      one_factorial[u] = Real(1) / std::tgamma(u + 1);
    return one_factorial;
  }();

  auto lambda_norm = [] {
    ffd::array1d<Real, order + 1> ret0;
    ret0.fill(0.);
    // std::cerr << "starting lambda_norm values:" << std::endl;
    ret0[0] = std::pow(std::abs(U_norm), order);
    for (ulong j = 1; j <= order; ++j) {
      ret0[j] = std::pow(std::abs(U_norm), j);
      // std::cerr << j << " " << ret0[j] << std::endl;
    }
    // ret0.fill(1.);
    return ret0;
  }();

  auto chosen_sets = [] {
    // ffd::array1d<ulong, actual_order + 1> ret;
    ffd::array1d<ulong, order + 1> ret;
    ret.fill(0);
    return ret;
  }();

  ffd::array2d<Real, order, order> M0;
  auto fill_G0_matrix_r = [&G0, &X_i, &M0]() mutable {
    auto psi_u = psi_f('u');
    auto bpsi_u = Bar(psi_f('u'));
    for (uint j = 0; j < order; ++j) {
      for (uint k = 0; k < order; ++k) {
        if (j != k) {
          M0(k, j) = -G0(std::make_pair(psi_u, X_i[j]),
                         std::make_pair(bpsi_u, X_i[k]));
        } else {
          M0(k, j) = 0;
        }
      }
    }
    // ffd::math_tools::print_matrix(M0);
    // exit(1);
  };

  ffd::vector3d<Complex> Sp(size_kw, actual_order + 1, (1ul << order));
  ffd::vector3d<Complex> Ch(size_kw, actual_order + 1, (1ul << order));
  ffd::array2d<Real, actual_order + 1, (1ul << order)> Density;
  ffd::array1d<Real, (1ul << order)> Weights;
  auto compute_weights_r = [&timers, &lambda_norm, size_k, size_w, size_kw,
                            &p_max, &M0, &exp_K, &Sp, &Ch, &Density,
                            &Weights]() mutable {
    ffd::array1d<poly_t, (1ul << actual_order)> Z, Sup, Sdn, Nup, Ndn, SupSdn,
        SupSup;
    ffd::array1d<poly_t, (1ul << order)> Z_, Nup_, Ndn_;
    // exit(1);

    ffd::vector1d<Complex> exp_kw(size_kw, 0.);

    Density.fill(0.);
    Sp.fill(0.);
    Ch.fill(0.);
    Weights.fill(0.);

    ffd::get<'p'>(timers).ini();
    auto PM = ffd::principal_minors::PrincipalMinorsBFS(M0);
    // for (BinaryInt S = 0; S < size(PM); ++S) {
    //   std::cout << S % (1 << actual_order) << " = "
    //             << std::abs(PM[S]) * (std::abs(PM[S]) > 1e-10) << std::endl;
    // }
    // exit(1);
    ffd::get<'p'>(timers).fin();

    ffd::get<'a'>(timers).ini();
    auto Zup = ffd::alpha_f::AlphaFunctionSet<order, actual_order>(PM);
    // auto Zup_ = ffd::alpha_f::AlphaFunction<order, true>(PM);
    auto Zdn = Zup;
    ffd::get<'a'>(timers).fin();
    // for (BinaryInt S = 0; S < size(Zup); ++S) {
    //   std::cout << S << "  Z = " << Zup[S][0] << "  " << Zup_[S][0]
    //             << std::endl;
    // }

    for (uint v1 = 0; v1 < order; ++v1) {
      const auto v1_ = (1 << v1);

      ffd::get<'a'>(timers).ini();
      const auto Aup =
          ffd::alpha_f::AlphaFunctionSet<order, actual_order>(PM, v1_);
      ffd::get<'a'>(timers).fin();

      const auto W = (1 << order) - 1 - v1_;
      const auto W_ = W + (1 << (order + 1));
      uint D = 0;
      for (uint D_ = W_; D_ != W; D_ = ((D_ - 1) & W_)) {
        const auto D__ = W_ - D_;
        const auto pop_D = ffd::core_math::popcount(D__);
        if (pop_D <= actual_order) {
          // Z[D] = Zup[D__] * Zdn[D__];
          // Nup[D] = Aup[D__] * Zdn[D__];
          Z_[D] = Zup[D__] * Zdn[D__];
          Nup_[D] = Aup[D__] * Zdn[D__];
          // if (D >= size(Z_)) {
          //   std::cout << "alert D! " << D << " " << size(Z) << std::endl;
          // }
          ffd::get<'r'>(timers).ini();
#pragma GCC unroll 4
#pragma GCC ivdep
          for (uint C = D; C != 0; C = ((C - 1) & D)) {
            // Nup[D] -= Nup[D - C] * Z[C];
            Nup_[D] -= Nup_[D - C] * Z_[C];
          }
          ffd::get<'r'>(timers).fin();
          for (uint r = 0; r <= p_max[pop_D]; ++r) {
            // Density(r, D__ + v1_) += std::pow(-1, pop_D) * 2. * Nup[D][r];
            Density(r, D__ + v1_) += std::pow(-1, pop_D) * 2. * Nup_[D][r];
          }
        }
        D++;
      }

      // for (BinaryInt S = 0; S < size(Zup); ++S) {
      //   std::cout << S << "  Z = " << Zup[S][0] << "  " << Zup_[S][0]
      //             << std::endl;
      // }
      // exit(1);

      for (uint v2 = 0; v2 < v1; ++v2) {
        const auto v2_ = (1 << v2);
        const auto v12_ = v1_ + v2_;
        const auto V = (1 << order) - 1 - v12_;
        const auto V_ = V + (1 << (order + 1));

        ffd::get<'a'>(timers).ini();
        const auto Adn =
            ffd::alpha_f::AlphaFunctionSet<order, actual_order>(PM, v2_);
        const auto Aupup =
            ffd::alpha_f::AlphaFunctionSet<order, actual_order>(PM, v12_);

        //   if (v1_ == 32 && v2_ == 16 && v12_ == 48) {
        //     std::cout << v1_ << " " << v2_ << " " << v12_ << std::endl;

        //     for (BinaryInt S = 0; S < size(Zup); ++S) {
        //       std::cout << S << "  Z = " << Zup[S][0] << std::endl;
        //     }
        //     for (BinaryInt S = 0; S < (1 << actual_order); ++S) {
        //       std::cout << S << " Aup = " << Aup[S][0] << std::endl;
        //     }
        //     for (BinaryInt S = 0; S < (1 << actual_order); ++S) {
        //       std::cout << S << " Adn = " << Adn[S][0] << std::endl;
        //     }
        //     for (BinaryInt S = 0; S < (1 << actual_order); ++S) {
        //       std::cout << S << " Aupup = " << Aupup[S][0] << std::endl;
        //     }
        //     exit(1);
        // }

        ffd::get<'a'>(timers).fin();

        for (uint w = 0; w < size_w; ++w) {
          for (uint k = 0; k < size_k; ++k) {
            exp_kw[k + w * size_k] = exp_K(v1, v2, w, k);
          }
        }

        uint S = 0;
        for (uint S_ = V_; S_ != V; S_ = ((S_ - 1) & V_)) {
          const auto S__ = V_ - S_;
          const auto pop_S = ffd::core_math::popcount(S__);

          Z[S] = Zup[S__] * Zdn[S__];
          Sdn[S] = Adn[S__] * Zup[S__];
          Nup[S] = Aup[S__] * Zdn[S__];
          SupSdn[S] = Aup[S__] * Adn[S__];
          SupSup[S] = Aupup[S__] * Zdn[S__];
          // if (S >= size(Z)) {
          //   std::cout << "alert S! " << S << " " << size(Z) << std::endl;
          // }
          ffd::get<'r'>(timers).ini();
#pragma GCC unroll 4
#pragma GCC ivdep
          for (uint C = S; C != 0; C = ((C - 1) & S)) {
            const auto S_C = S - C;
            Nup[S] -= Nup[S_C] * Z[C];
            const auto NupZSdnZ = Nup[S_C] * Sdn[C];
            // SupSdn[S] -= SupSdn[S_C] * Z[C];
            // SupSup[S] -= SupSup[S_C] * Z[C];
            SupSdn[S] -= SupSdn[S_C] * Z[C] + NupZSdnZ;
            SupSup[S] -= SupSup[S_C] * Z[C] + NupZSdnZ;
          }
          SupSdn[S] -= Nup[0] * Sdn[S];
          SupSup[S] -= Nup[0] * Sdn[S];
          ffd::get<'r'>(timers).fin();

          // const auto Sp_ = SupSup[S];
          // const auto Ch_ = SupSdn[S];
          const auto Sp_ = (SupSup[S] - SupSdn[S]) * 2.;
          const auto Ch_ = (SupSup[S] + SupSdn[S]) * 2.;

          ffd::get<'w'>(timers).ini();
          Real pow_ashift = 1.;
          for (uint r = 0; r <= p_max[pop_S]; ++r) {
            // Weights[S__ + v12_] += std::abs(Sp_[r]) * pow_ashift;
            Weights[S__ + v12_] +=
                (std::abs(Sp_[r]) + std::abs(Ch_[r])) * pow_ashift;
            const auto pos = size_kw * (r + (S__ + v12_) * (actual_order + 1));
#pragma GCC unroll 4
#pragma GCC ivdep
            for (uint kw = 0; kw < size_kw; ++kw) {
              Sp[kw + pos] += exp_kw[kw] * Sp_[r];
              Ch[kw + pos] += exp_kw[kw] * Ch_[r];
            }
            pow_ashift *= ashift_target;
          }
          ffd::get<'w'>(timers).fin();
          S++;
        }
        // exit(1);
      }
    }
    Weights[0] = lambda_norm[0];
    for (uint S = 3; S < (1 << order); ++S) {
      const auto pop_S = ffd::core_math::popcount(S);
      if (pop_S >= 2) {
        const auto pop_S_1 = pop_S - 1;
        const auto pop_S_2 = pop_S - 2;
        const auto one_factor = 1. / (0.5 * pop_S * (pop_S - 1));
        const auto den_factor = (pop_S != order) ? 1. / (pop_S) : 1.;
        for (uint r = 0; r <= p_max[pop_S_1]; ++r) {
          Density(r, S) *= den_factor;
        }
        Weights[S] *= lambda_norm[pop_S] * one_factor;
        for (uint r = 0; r <= p_max[pop_S_2]; ++r) {
          const auto pos = size_kw * (r + S * (actual_order + 1));
          for (uint kw = 0; kw < size_kw; ++kw) {
            Sp[kw + pos] *= one_factor;
            Ch[kw + pos] *= one_factor;
          }
        }
      }
    }
    // exit(1);
  };

  auto [Sp_fnames, Ch_fnames, Den_fname] = [k_vec, thread_ID = thread_ID,
                                            k_bravais = k_bravais, mu0 = mu0] {
    std::stringstream ss;
    if (geometry == square)
      ss << "square_";
    if (geometry == triangle)
      ss << "triang_";
    ss << Lx << "x" << Ly << "_B" << Beta << "_n0-" << n0 << "_t" << t_hopping
       << "_tp" << tprime_hopping;
    fs::create_directory(ss.str());
    fs::current_path(ss.str());
    fs::remove("kill");
    fs::create_directory("parameters");
    fs::current_path("parameters");
    {
      std::ofstream Lx_file("Lx");
      Lx_file << Lx;
    }
    {
      std::ofstream Ly_file("Ly");
      Ly_file << Ly;
    }
    {
      std::ofstream geometry_file("geometry");
      geometry_file << geometry;
    }
    {
      std::ofstream Beta_file("Beta");
      Beta_file << Beta;
    }
    {
      std::ofstream n0_file("n0");
      n0_file << n0;
    }
    {
      std::ofstream t_hopping_file("t_hopping");
      t_hopping_file << t_hopping;
    }
    {
      std::ofstream tprime_hopping_file("tprime_hopping");
      tprime_hopping_file << tprime_hopping;
    }
    {
      std::ofstream mu0_file("mu0");
      mu0_file << std::setprecision(15) << mu0;
    }
    fs::current_path("..");
    fs::create_directory("density");
    std::stringstream ss_density;
    ss_density << "density/" << thread_ID;
    ffd::vector3d<std::string> ret_Sp(2, omega_end - omega_begin, size(k_vec));
    ffd::vector3d<std::string> ret_Ch(2, omega_end - omega_begin, size(k_vec));
    for (ulong k = 0; k < size(k_vec); ++k) {
      for (ulong w = omega_begin; w < omega_end; ++w) {
        std::stringstream ss_Sp_i, ss_Sp_r, ss_Ch_i, ss_Ch_r, ss_Sp, ss_Ch,
            ss_Sp_r2, ss_Ch_r2, ss_Sp_i2, ss_Ch_i2;
        auto [kx, ky] = k_vec[k];
        Real Kx = kx * k_bravais[0][0] + ky * k_bravais[1][0];
        if (std::abs(Kx) < 1e-13)
          Kx = 0;
        Real Ky = kx * k_bravais[0][1] + ky * k_bravais[1][1];
        if (std::abs(Ky) < 1e-13)
          Ky = 0;
        ss_Sp << "Sp_"
              << "kx" << Kx << "_ky" << Ky << "_om" << w;
        ss_Ch << "Ch_"
              << "kx" << Kx << "_ky" << Ky << "_om" << w;
        ss_Sp_r << ss_Sp.str() << "_r";
        ss_Sp_i << ss_Sp.str() << "_i";
        ss_Ch_r << ss_Ch.str() << "_r";
        ss_Ch_i << ss_Ch.str() << "_i";
        fs::create_directory(ss_Sp_r.str());
        fs::create_directory(ss_Sp_i.str());
        fs::create_directory(ss_Ch_r.str());
        fs::create_directory(ss_Ch_i.str());
        {
          ss_Sp_r2 << ss_Sp_r.str();
          ss_Sp_r2 << "/parameters";
          fs::create_directory(ss_Sp_r2.str());
          fs::current_path(ss_Sp_r2.str());
          {
            std::ofstream kx_file("kx");
            kx_file << std::setprecision(15) << Kx;
          }
          {
            std::ofstream ky_file("ky");
            ky_file << std::setprecision(15) << Ky;
          }
          {
            std::ofstream omega_file("omega");
            omega_file << w;
          }
          {
            std::ofstream component_file("component");
            component_file << 'r';
          }
          fs::current_path("..");
          fs::current_path("..");
        }
        {
          ss_Ch_r2 << ss_Ch_r.str();
          ss_Ch_r2 << "/parameters";
          fs::create_directory(ss_Ch_r2.str());
          fs::current_path(ss_Ch_r2.str());
          {
            std::ofstream kx_file("kx");
            kx_file << std::setprecision(15) << Kx;
          }
          {
            std::ofstream ky_file("ky");
            ky_file << std::setprecision(15) << Ky;
          }
          {
            std::ofstream omega_file("omega");
            omega_file << w;
          }
          {
            std::ofstream component_file("component");
            component_file << 'r';
          }
          fs::current_path("..");
          fs::current_path("..");
        }
        {
          ss_Sp_i2 << ss_Sp_i.str();
          ss_Sp_i2 << "/parameters";
          fs::create_directory(ss_Sp_i2.str());
          fs::current_path(ss_Sp_i2.str());
          {
            std::ofstream kx_file("kx");
            kx_file << std::setprecision(15) << Kx;
          }
          {
            std::ofstream ky_file("ky");
            ky_file << std::setprecision(15) << Ky;
          }
          {
            std::ofstream omega_file("omega");
            omega_file << w;
          }
          {
            std::ofstream component_file("component");
            component_file << 'i';
          }
          fs::current_path("..");
          fs::current_path("..");
        }
        {
          ss_Ch_i2 << ss_Ch_i.str();
          ss_Ch_i2 << "/parameters";
          fs::create_directory(ss_Ch_i2.str());
          fs::current_path(ss_Ch_i2.str());
          {
            std::ofstream kx_file("kx");
            kx_file << std::setprecision(15) << Kx;
          }
          {
            std::ofstream ky_file("ky");
            ky_file << std::setprecision(15) << Ky;
          }
          {
            std::ofstream omega_file("omega");
            omega_file << w;
          }
          {
            std::ofstream component_file("component");
            component_file << 'i';
          }
          fs::current_path("..");
          fs::current_path("..");
        }
        // fs::create_directory(ss_r2.str());
        ss_Sp_r << "/" << thread_ID;
        ss_Sp_i << "/" << thread_ID;
        ss_Ch_r << "/" << thread_ID;
        ss_Ch_i << "/" << thread_ID;
        ret_Sp(0, w, k) = ss_Sp_r.str();
        ret_Sp(1, w, k) = ss_Sp_i.str();
        ret_Ch(0, w, k) = ss_Ch_r.str();
        ret_Ch(1, w, k) = ss_Ch_i.str();
      }
    }
    return std::make_tuple(ret_Sp, ret_Ch, ss_density.str());
  }();

  auto set_lambda_norm_r = [&chosen_sets, &lambda_norm, &normalization_phase,
                            &timers, tid = thread_ID,
                            pid = process_ID]() mutable {
    static_assert(actual_order > 2);
    auto const t_phase_1 = (0.7 + 0.7 * Proba()) * thermalization_time;
    auto const t_phase_2 = (3. + Proba()) * thermalization_time;

    if (normalization_phase != 0) {
      if (ffd::get<'T'>(timers)() > t_phase_1 && normalization_phase == 1) {
        normalization_phase = 2;
        {
          std::stringstream file_ss;
          for (ulong j = 0; j < order; ++j) {
            file_ss << lambda_norm[j] << "\n";
          }
          if (verbose_stdout >= 1)
            std::cout << file_ss.str();
          std::ofstream file_lambda("../" + pid + "/lambda_norm/" + tid);
          file_lambda << file_ss.str();
        }
      } else if (ffd::get<'T'>(timers)() > t_phase_2 &&
                 normalization_phase == 2) {
        normalization_phase = 0;
        std::array<std::vector<Real>, order> lambda_norm_file;
        for (auto& p : fs::directory_iterator("../" + pid + "/lambda_norm")) {
          std::ifstream ff(p.path());
          for (ulong j = 0; j < order; ++j) {
            Real val;
            ff >> val;
            lambda_norm_file[j].push_back(val);
          }
        }
        for (ulong j = 0; j < order; ++j) {
          std::nth_element(
              lambda_norm_file[j].begin(),
              lambda_norm_file[j].begin() + lambda_norm_file[j].size() / 2,
              lambda_norm_file[j].end());
          lambda_norm[j] = lambda_norm_file[j][lambda_norm_file[j].size() / 2];
          if (verbose_stdout >= 1)
            std::cerr << j << " lambda_norm_median = " << lambda_norm[j]
                      << "\n";
        }
      } else {
        auto const [opt_times, act_times] = [&chosen_sets] {
          ffd::array1d<Real, order> opt_times;
          opt_times.fill((1. - time_normalization_sector - time_target_order) /
                         (order - 2));
          opt_times[0] = time_normalization_sector;
          opt_times[order - 1] = time_target_order;
          ffd::array1d<Real, order> act_times;
          Real acc = 0;
          for (ulong j = 0; j < order; ++j) {
            act_times[j] = chosen_sets[j] + 1;
            acc += act_times[j];
          }
          for (ulong j = 0; j < order; ++j) {
            act_times[j] /= acc;
          }
          return std::make_pair(opt_times, act_times);
        }();
        for (ulong j = 0; j < order; ++j) {
          if (j == 0 || j == order - 1) {
            lambda_norm[j] *=
                std::exp(thermalization_speed *
                         std::log(opt_times[j] / act_times[j]) / order);
          }
        }
      }
    }
  };

  accum_t Norm;
  ffd::array2d<accum_t, actual_order + 2, actual_order + 2> Den_acc;
  ffd::vector3d<accum_t> Sp_acc_re(size_kw, actual_order + 1, actual_order + 1);
  ffd::vector3d<accum_t> Sp_acc_im(size_kw, actual_order + 1, actual_order + 1);
  ffd::vector3d<accum_t> Ch_acc_re(size_kw, actual_order + 1, actual_order + 1);
  ffd::vector3d<accum_t> Ch_acc_im(size_kw, actual_order + 1, actual_order + 1);
  auto accumulate_r = [size_k, size_w, size_kw, &p_max, &normalization_phase,
                       &Weights, &Density, &Den_acc, &Sp, &Ch, &Sp_acc_re,
                       &Sp_acc_im, &Ch_acc_re, &Ch_acc_im,
                       &Norm](auto const& times, auto const t_MC) mutable {
    if (normalization_phase != 0)
      return;

    Norm += times[0] / Weights[0];
    for (uint S = 1; S < (1 << order); ++S) {
      const auto pop_S = ffd::core_math::popcount(S);
      const auto pop_S_1 = pop_S - 1;
      const auto pop_S_2 = pop_S - 2;
      const auto one_Weight_S =
          std::abs(Weights[S]) != 0 ? times[S] / Weights[S] : Real(0);
      for (ulong r = 0; r <= p_max[pop_S_1]; ++r) {
        Den_acc(r, pop_S_1) &= Density(r, S) * one_Weight_S;
      }
      if (pop_S > 1) {
        for (ulong r = 0; r <= p_max[pop_S_2]; ++r) {
          for (ulong kw = 0; kw < size_kw; ++kw) {
            Sp_acc_re(kw, r, pop_S_2) &= real(Sp(kw, r, S)) * one_Weight_S;
            Sp_acc_im(kw, r, pop_S_2) &= imag(Sp(kw, r, S)) * one_Weight_S;
            Ch_acc_re(kw, r, pop_S_2) &= real(Ch(kw, r, S)) * one_Weight_S;
            Ch_acc_im(kw, r, pop_S_2) &= imag(Ch(kw, r, S)) * one_Weight_S;
          }
        }
      }
    }
    for (uint j = 0; j <= actual_order; ++j) {
      for (ulong r = 0; r <= p_max[j]; ++r) {
        Den_acc(r, j).exec_deferred_add();
        for (ulong kw = 0; kw < size_kw; ++kw) {
          Sp_acc_re(kw, r, j).exec_deferred_add();
          Sp_acc_im(kw, r, j).exec_deferred_add();
          Ch_acc_re(kw, r, j).exec_deferred_add();
          Ch_acc_im(kw, r, j).exec_deferred_add();
        }
      }
    }
    // exit(1);
  };

  std::string log_file_name{"../" + process_ID + "/log.log"};
  auto print_timers_r = [&timers, &chosen_sets, log_file_name,
                         t_counter = (1. + Proba()) * print_interval,
                         &lambda_norm](ulong t_MC = 0) mutable {
    if (ffd::get<'T'>(timers)() > t_counter) {
      std::stringstream ss;
      auto print_timer = [&ss, &timers, t_MC](auto t,
                                              std::string const& x) mutable {
        ss << "t(" << x << ")= " << 1000 * 1000 * t() / t_MC << " \u03BCs   "
           << 100 * t() / ffd::get<'T'>(timers)() << " %\n";
        t.reset();
      };
      print_timer(ffd::get<'M'>(timers), "matrix  ");
      print_timer(ffd::get<'E'>(timers), "exp     ");
      print_timer(ffd::get<'p'>(timers), "minors  ");
      print_timer(ffd::get<'a'>(timers), "alpha   ");
      print_timer(ffd::get<'r'>(timers), "recurs  ");
      print_timer(ffd::get<'w'>(timers), "k-omega ");
      print_timer(ffd::get<'H'>(timers), "heatbath");
      print_timer(ffd::get<'K'>(timers), "accum   ");
      print_timer(ffd::get<'N'>(timers), "norm    ");
      print_timer(ffd::get<'T'>(timers), "total   ");
      for (ulong j = 0; j <= order; ++j)
        ss << "card=" << j << ", t=" << chosen_sets[j] * 100. / t_MC
           << "%, lambda_norm=" << lambda_norm[j] << "\n";
      if (t_MC != 0)
        ss << "after " << t_MC << " iterations "
           << "\n";
      if (verbose_stdout >= 2)
        std::cout << ss.str();
      if (verbose_file >= 2) {
        std::ofstream out_log(log_file_name);
        out_log << ss.str();
      }
      t_counter += (1 + Proba()) * print_interval;
    }
  };

  auto print_accumulators_r =
      [size_k, size_w, size_kw, p_max, &normalization_phase, &timers,
       t_counter = (1. + Proba()) * print_interval, &Sp_acc_re, &Sp_acc_im,
       &Ch_acc_re, &Ch_acc_im, &Den_acc, log_file_name, one_factorial, &Norm,
       Sp_fnames = Sp_fnames, Ch_fnames = Ch_fnames,
       Den_fname = Den_fname]() mutable {
        using std::setprecision;
        if (ffd::get<'T'>(timers)() > t_counter && !normalization_phase) {
          auto const norm_val = ENumber(Norm.mean_error());

          std::stringstream to_stdout, to_density_file, first_line_file,
              Sp_to_stdout, Ch_to_stdout, W_to_stdout, Norm_to_stdout;

          Norm_to_stdout << "Norm= " << setprecision(10) << to_string(norm_val)
                         << "\n";
          first_line_file << -1 << " " << -1 << " " << setprecision(10)
                          << norm_val.value << " " << setprecision(3)
                          << norm_val.error << "\n";

          to_density_file << 0 << " " << 0 << " " << n0 << " " << 0 << "\n";
          for (ulong j = 1; j <= actual_order; ++j) {
            for (ulong p = 0; p <= p_max[j]; ++p) {
              auto n_err = ENumber(Den_acc(p, j).mean_error());
              n_err *= one_factorial[j];
              auto const n_err_norm = n_err / norm_val;
              to_stdout << "n_{" << j << ", " << p
                        << "} = " << to_string(n_err_norm) << "\n";
              to_density_file << j << " " << p << " " << setprecision(8)
                              << n_err.value << " " << setprecision(3)
                              << n_err.error << "\n";
            }
          }
          {
            std::ofstream density_f(Den_fname, std::ofstream::trunc);
            density_f << first_line_file.str() << to_density_file.str();
          }
          for (ulong k = 0; k < size_k; ++k) {
            for (ulong w = 0; w < size_w; ++w) {
              auto kw = k + w * size_k;
              std::stringstream Sp_to_file_r, Sp_to_file_i, Ch_to_file_r,
                  Ch_to_file_i;
              for (ulong j = 0; j <= actual_order; ++j) {
                for (ulong p = 0; p <= p_max[j]; ++p) {
                  auto Sp_r = ENumber(Sp_acc_re(kw, p, j).mean_error());
                  auto Sp_i = ENumber(Sp_acc_im(kw, p, j).mean_error());
                  auto Ch_r = ENumber(Ch_acc_re(kw, p, j).mean_error());
                  auto Ch_i = ENumber(Ch_acc_im(kw, p, j).mean_error());
                  Sp_r *= one_factorial[j] * pow(-1, j);
                  Sp_i *= one_factorial[j] * pow(-1, j);
                  Ch_r *= one_factorial[j] * pow(-1, j);
                  Ch_i *= one_factorial[j] * pow(-1, j);
                  auto const Sp_r_norm = Sp_r / norm_val;
                  auto const Sp_i_norm = Sp_i / norm_val;
                  auto const Ch_r_norm = Ch_r / norm_val;
                  auto const Ch_i_norm = Ch_i / norm_val;
                  Sp_to_stdout << Sp_fnames(0, w, k);
                  Ch_to_stdout << Ch_fnames(0, w, k);
                  Sp_to_stdout << ", " << j << " " << p << " = "
                               << to_string(Sp_r_norm) << " +I*( "
                               << to_string(Sp_i_norm) << " )\n";
                  Ch_to_stdout << ", " << j << " " << p << " = "
                               << to_string(Ch_r_norm) << " +I*( "
                               << to_string(Ch_i_norm) << " )\n";
                  Sp_to_file_r << j << " " << p << " " << setprecision(8)
                               << Sp_r.value << " " << setprecision(3)
                               << Sp_r.error << "\n";
                  Sp_to_file_i << j << " " << p << " " << setprecision(8)
                               << Sp_i.value << " " << setprecision(3)
                               << Sp_i.error << "\n";
                  Ch_to_file_r << j << " " << p << " " << setprecision(8)
                               << Ch_r.value << " " << setprecision(3)
                               << Ch_r.error << "\n";
                  Ch_to_file_i << j << " " << p << " " << setprecision(8)
                               << Ch_i.value << " " << setprecision(3)
                               << Ch_i.error << "\n";
                }
              }
              {
                std::stringstream ss;
                ss << first_line_file.str() << Sp_to_file_r.str();
                std::ofstream Sp_f_r(Sp_fnames(0, w, k), std::ofstream::trunc);
                Sp_f_r << ss.str();
              }
              {
                std::stringstream ss;
                ss << first_line_file.str() << Sp_to_file_i.str();
                std::ofstream Sp_f_i(Sp_fnames(1, w, k), std::ofstream::trunc);
                Sp_f_i << ss.str();
              }
              {
                std::stringstream ss;
                ss << first_line_file.str() << Ch_to_file_r.str();
                std::ofstream Ch_f_r(Ch_fnames(0, w, k), std::ofstream::trunc);
                Ch_f_r << ss.str();
              }
              {
                std::stringstream ss;
                ss << first_line_file.str() << Ch_to_file_i.str();
                std::ofstream Ch_f_i(Ch_fnames(1, w, k), std::ofstream::trunc);
                Ch_f_i << ss.str();
              }
            }
          }
          if (verbose_stdout >= 3) {
            std::cout << Norm_to_stdout.str();
            std::cout << to_stdout.str();
            std::cout << Sp_to_stdout.str();
            std::cout << Ch_to_stdout.str();
          }
          t_counter += (1 + Proba()) * print_interval;
        }
      };

  auto check_for_kill_signal_r = [t_counter =
                                      (1. + Proba()) * kill_signal_interval,
                                  &timers, process_ID = process_ID]() mutable {
    if (ffd::get<'T'>(timers)() > t_counter) {
      if (fs::exists("kill") || fs::exists("../kill") ||
          fs::exists(process_ID + "/kill") || fs::exists("../../kill"))
        std::exit(0);
      t_counter += (1 + Proba()) * kill_signal_interval;
    }
  };

  // MONTE CARLO STARTS HERE
  for (ulong t_MC = 0;; ++t_MC) {
    ffd::get<'T'>(timers).ini();
    ffd::get<'M'>(timers).ini();
    fill_G0_matrix_r();
    ffd::get<'M'>(timers).fin();
    ffd::get<'E'>(timers).ini();
    exp_K_r();
    ffd::get<'E'>(timers).fin();
    ffd::get<'W'>(timers).ini();
    compute_weights_r();
    ffd::get<'W'>(timers).fin();
    ffd::get<'H'>(timers).ini();
    auto const [times, X_new, chosen_set] =
        ffd::heat_bath_mc::HeatBathAll<order>(Weights, proposer, X_i);
    ffd::get<'H'>(timers).fin();
    ffd::get<'N'>(timers).ini();
    ++chosen_sets[ffd::core_math::popcount(chosen_set)];
    set_lambda_norm_r();
    ffd::get<'N'>(timers).fin();
    ffd::get<'K'>(timers).ini();
    accumulate_r(times, t_MC);
    ffd::get<'K'>(timers).fin();
    X_i = X_new;
    ffd::get<'T'>(timers).fin();
    print_timers_r(t_MC + 1);
    print_accumulators_r();
    check_for_kill_signal_r();
  }
}  // main