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

  auto const G0 = [H0 = H0] {
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

  // std::cerr << "k_bravais0 = " << std::setprecision(16) << k_bravais[0][0]
  //           << " " << k_bravais[0][1] << "\n";
  // std::cerr << "k_bravais1 = " << k_bravais[1][0] << " " << k_bravais[1][1]
  //           << "\n";

  auto X = Spacetime_g<d>();
  auto X_i = [X] {
    std::array<std::decay_t<decltype(X(0., 0, 0))>, order> X_i;
    for (uint j = 0; j < order; ++j)
      X_i[j] = X(Beta * Proba(), 0, 0);
    return X_i;
  }();

  std::size_t constexpr poly_n = order;  // (order + 1) / 2 - 1;
  using poly_t = ffd::taylor_polynomial::P<Real, poly_n>;
  using mat_el_t = poly_t;  // Real;  // poly_t;
  auto const alpha = [] {
    poly_t alpha{0};
    alpha[1] = Real(1);
    if constexpr (std::is_same_v<mat_el_t, poly_t>) {
      return alpha;
    } else {
      return Real(0);
    }
  }();
  ffd::array2d<mat_el_t, order, order> M;
  ffd::array2d<Real, order, order> M0;
  auto fill_G0_matrix_r = [psi_u = psi_f('u'), bpsi_u = Bar(psi_f('u')), alpha,
                           &G0, &X_i, &M, &M0]() mutable {
    //    M.fill(Real(0));
    M0.fill(Real(0));
    for (uint j = 0; j < size_y(M); ++j) {
      for (uint k = 0; k < size_x(M); ++k) {
        if (j != k) {
          M0(k, j) = -G0(std::make_pair(psi_u, X_i[j]),
                         std::make_pair(bpsi_u, X_i[k]));
          //          M(k, j) = M0(k, j);
        } else {
          //          M(j, k) = alpha;
          M0(k, j) = 0;
        }
      }
    }
  };

  auto lattice_equivalent = [](std::array<int, 2> n0,
                               std::array<int, 2> n1) -> bool {
    static_assert(Lx == Ly);
    auto constexpr L = Lx;
    auto are_equivalent = [](auto n2, auto n3) -> bool {
      bool ret = true;
      for (ulong s = 0; s < 2; ++s) {
        int diff = n2[s] - n3[s];
        while (diff >= L)
          diff -= L;
        while (diff < 0)
          diff += L;
        ret = ret && diff == 0;
      }  // for s in range(0, 2)
      return ret;
    };
    auto swap_k = [](auto& n) {
      auto t = n[0];
      n[0] = n[1];
      n[1] = t;
    };
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
    return false;
  };

  auto const k_vec = [lattice_equivalent] {
    static_assert(Lx == Ly);
    std::vector<std::array<int, 2>> vector_k;
    for (int fac = 1; fac <= Lx; ++fac) {
      if (Lx % fac == 0) {
        vector_k.resize(0);
        for (int ny = 0; ny < Ly; ny += fac) {
          for (int nx = 0; nx < Lx; nx += fac) {
            std::array<int, 2> candidate{nx, ny};
            bool not_inside = true;
            for (auto v : vector_k)
              not_inside = not_inside && !lattice_equivalent(candidate, v);

            bool pb = false;
            pb |= (k_slice == QQ && candidate[0] == candidate[1]);
            pb |= (k_slice == QP && candidate[0] == Lx / 2);
            pb |= (k_slice == Q0 && candidate[1] == 0.);
            pb |= (k_slice == BZ);
            if (not_inside && pb)
              vector_k.push_back(candidate);
          }  // for nx in range(0, size_x(ret))
        }    // for ny in range(0, size_y(ret))
        if (size(vector_k) <= n_k_points_max)
          break;
      }
    }  // for fac in range(1, Lx)
    std::vector<std::array<Real, 2>> ret;
    for (ulong j = 0; j < size(vector_k); ++j) {
      ret.push_back(
          {2 * vector_k[j][0] * M_PI / Lx, 2 * vector_k[j][1] * M_PI / Ly});
      // std::cout << j << "  " << vector_k[j][0] << " " << vector_k[j][1]
      //           << std::endl;
    }  // for j in range(0, size(vector_k))
    if (verbose_stdout >= 1)
      std::cout << "size_k_vec = " << size(ret) << "\n";
    return ret;
  }();
  //  exit(1);

  auto [Sigma_fnames, density_fname] = [k_vec, thread_ID = thread_ID,
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
    {
      ffd::external::json20::dictionary jpar;
      jpar["Lx"] = Lx;
      jpar["Ly"] = Ly;
      jpar["geometry"] = geometry;
      jpar["Beta"] = Beta;
      jpar["n0"] = n0;
      jpar["t_hopping"] = t_hopping;
      jpar["tprime_hopping"] = tprime_hopping;
      jpar["mu0"] = mu0;
      ffd::external::json20::to_file("parameters.json", jpar);
    }
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
    ffd::vector3d<std::string> ret0(2, omega_end - omega_begin, size(k_vec));
    for (ulong k = 0; k < size(k_vec); ++k) {
      for (ulong j = 0; j < omega_end - omega_begin; ++j) {
        std::stringstream ss_i, ss_r, ss_f, ss_r2, ss_i2;
        auto [kx, ky] = k_vec[k];
        Real Kx = kx * k_bravais[0][0] + ky * k_bravais[1][0];
        if (std::abs(Kx) < 1e-13)
          Kx = 0;
        Real Ky = kx * k_bravais[0][1] + ky * k_bravais[1][1];
        if (std::abs(Ky) < 1e-13)
          Ky = 0;
        ss_f << "kx" << Kx << "_ky" << Ky << "_om" << j + omega_begin;
        ss_r << ss_f.str() << "_r";
        ss_i << ss_f.str() << "_i";
        fs::create_directory(ss_r.str());
        fs::create_directory(ss_i.str());
        ss_r2 << ss_r.str();
        ss_r2 << "/parameters";
        fs::create_directory(ss_r2.str());
        fs::current_path(ss_r2.str());
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
          omega_file << j + omega_begin;
        }
        {
          std::ofstream component_file("component");
          component_file << 'r';
        }
        fs::current_path("..");
        fs::current_path("..");
        ss_i2 << ss_i.str();
        ss_i2 << "/parameters";
        fs::create_directory(ss_i2.str());
        fs::current_path(ss_i2.str());
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
          omega_file << j + omega_begin;
        }
        {
          std::ofstream component_file("component");
          component_file << 'i';
        }
        fs::current_path("..");
        fs::current_path("..");
        fs::create_directory(ss_r2.str());
        ss_r << "/" << thread_ID;
        ss_i << "/" << thread_ID;
        ret0(0, j, k) = ss_r.str();
        ret0(1, j, k) = ss_i.str();
      }  // for k in range(0, size(vec))
    }    // for j in range(omega_begin, omega_end)
    return std::make_pair(ret0, ss_density.str());
  }();

  auto lambda_norm = [] {
    ffd::array1d<Real, order + 1> ret0;
    for (ulong j = 0; j < order + 1; ++j) {
      ret0[j] = std::pow(std::abs(U_norm), j + order * (j == 0));
    }  // for j in range(0, order+1)
    return ret0;
  }();
  int normalization_phase = 1;

  auto chosen_sets = [] {
    ffd::array1d<ulong, order + 1> ret;
    ret.fill(0);
    return ret;
  }();

  auto set_lambda_norm_r = [&chosen_sets, &lambda_norm, &normalization_phase,
                            &timers, tid = thread_ID, pid = process_ID,
                            t_phase_2 = (3 + Proba()) * thermalization_time,
                            t_phase_1 = (0.7 + 0.7 * Proba()) *
                                        thermalization_time]() mutable {
    static_assert(order > 2);
    if (normalization_phase != 0) {
      if (ffd::get<'T'>(timers)() > t_phase_1 && normalization_phase == 1) {
        normalization_phase = 2;
        {
          std::stringstream file_ss;
          for (ulong j = 0; j < order; ++j) {
            file_ss << lambda_norm[j] << "\n";
          }  // for j in range(0, order)
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
          }  // for j in range(0, order)
        }
        for (ulong j = 0; j < order; ++j) {
          std::nth_element(
              lambda_norm_file[j].begin(),
              lambda_norm_file[j].begin() + lambda_norm_file[j].size() / 2,
              lambda_norm_file[j].end());
          lambda_norm[j] = lambda_norm_file[j][lambda_norm_file[j].size() / 2];
          if (verbose_stdout >= 1)
            std::cout << j << " lambda_norm_median = " << lambda_norm[j]
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
          if (j == 0 || j == order - 1)
            lambda_norm[j] *=
                std::exp(thermalization_speed *
                         std::log(opt_times[j] / act_times[j]) / order);
        }
      }
    }
  };

  ffd::vector3d<poly_t> Sigma(order, order, (1ul << order));
  ffd::vector3d<poly_t> rho_T = Sigma;
  auto const bitmap = ffd::sigmadet_hubbard::CDet_Bitmap<order>();
  ffd::vector1d<Real> Sigma_a((1ul << order), 0);
  auto constexpr p_max_n = []() constexpr {
    std::array<ulong, order + 1> p_max = {0};
    for (int j = 0; j < order + 1; ++j)
      p_max[j] = std::min<int>(j - (j > 0), order - j);
    return p_max;
  }
  ();
  auto constexpr p_max_S = []() constexpr {
    std::array<ulong, order + 1> p_max = {0};
    for (int j = 0; j < order + 1; ++j)
      p_max[j] = std::min<int>(j - (j > 0) - (j > 1), order - j);
    return p_max;
  }
  ();

  auto Sigma_a_r = [p_max_n, p_max_S, &bitmap, &Sigma, &Sigma_a,
                    &lambda_norm]() mutable {
    using std::abs;
    Sigma_a(0) = lambda_norm(0);
    for (ulong S = 1; S < (1ul << order); ++S) {
      auto const bm_0_S = bitmap(0, S);
      Real ret = 0.;
      for (ulong j = 0; j < bm_0_S; ++j) {
        auto const bm_j = bitmap(j + 1, S);
        for (ulong l = 0; l < bm_0_S; ++l) {
          auto const bm_l = bitmap(l + 1, S);
          Real pow_ashift = 1;
          auto const p_max = bm_j == bm_l ? p_max_n[bm_0_S] : p_max_S[bm_0_S];
          for (ulong p = 0; p <= p_max; ++p, pow_ashift *= ashift_target) {
            ret += abs(Sigma(bm_l, bm_j, S)[p] * pow_ashift);
          }  // for p in range(0, poly_n+1)
        }    // for l in range(0, bm_0_S)
      }      // for j in range(0, bm_0_S)
      Sigma_a(S) = ret * lambda_norm(bm_0_S) / (bm_0_S * bm_0_S);
    }  // for S in range(0, (1ul<<order))
    return;
  };

  ffd::vector4d<Complex> exp_K(order, order, omega_end - omega_begin,
                               size(k_vec), Real(0));
  auto exp_K_r = [&normalization_phase, &k_vec, &exp_K, &X_i]() mutable {
    if (  // (t_MC%accumulate_every) != 0 ||
        normalization_phase != 0)
      return;
    ffd::array2d<Complex, order, omega_end - omega_begin> exp_om;
    for (ulong j = 0; j < size(exp_om); ++j) {
      auto const [u, om] = exp_om.indexes(j);
      exp_om[j] = std::exp(Complex(
          0., -X_i[u].first * (om + omega_begin + .5) * 2. * M_PI / Beta));
    }  // for j in sizerange(exp_om)

    auto exp_k = ffd::vector2d<Complex>{order, size(k_vec)};
    for (ulong j = 0; j < size(exp_k); ++j) {
      auto const [u, k] = exp_k.indexes(j);
      auto const [kx, ky] = k_vec[k];
      auto const [x, y] = X_i[u].second;
      exp_k[j] = std::exp(Complex(0., x * kx + y * ky));
    }  // for j in sizerange(exp_k)

    auto exp_om_k =
        ffd::vector3d<Complex>{order, omega_end - omega_begin, size(k_vec)};
    for (ulong j = 0; j < size(exp_om_k); ++j) {
      auto const [u, om, k] = exp_om_k.indexes(j);
      exp_om_k[j] = exp_k(u, k) * exp_om(u, om);
    }  // for j in sizerange(exp_om_k)

    for (ulong j = 0; j < size(exp_K); ++j) {
      auto const [u, v, om, k] = exp_K.indexes(j);
      if (u != v) {
        exp_K[j] = exp_om_k(u, om, k) * conj(exp_om_k(v, om, k));
      }
    }  // for j in sizerange(exp_K)
    return;
  };

  using accum_t = ffd::blocking_method::block_t<1024, 16, Real>;
  ffd::vector5d<accum_t> Sigma_acc(2, poly_n + 1, omega_end - omega_begin,
                                   size(k_vec), order + 1);
  ffd::array2d<accum_t, poly_n + 1, order + 1> Un_acc;
  auto Sigma_K_acc_r = [p_max_n, p_max_S, &normalization_phase, &bitmap,
                        &Un_acc, &Sigma, &exp_K, &Sigma_a, &Sigma_acc](
                           auto const& times, auto const t_MC) mutable {
    if (normalization_phase != 0)
      return;
    Sigma_acc(0, 0, 0, 0, 0) += times[0] / Sigma_a(0);
    for (ulong S = 1; S < (1ul << order); ++S) {
      auto const bm_0_S = bitmap(0, S);
      auto const one_Sigma_a_S =
          std::abs(Sigma_a[S]) != 0 ? times[S] / Sigma_a(S) : Real(0);
      for (ulong u = 0; u < bm_0_S; ++u) {
        auto const bm_u = bitmap(u + 1, S);
        for (ulong p = 0; p <= p_max_n[bm_0_S]; ++p) {
          Un_acc(p, bm_0_S) &= Sigma(bm_u, bm_u, S)[p] * one_Sigma_a_S;
        }  // for p in range(0, n_poly+1)
      }    // for u in range(0, order)
    }
    for (ulong j = 0; j < size(Un_acc); ++j) {
      Un_acc[j].exec_deferred_add();
    }  // for j in sizerange(Un_acc)
    auto constexpr order2 = order * order;
    {
      auto const V = (1ul << order) - 1;
      auto const index_V = V * order2;
      auto const one_Sigma_a_V =
          std::abs(Sigma_a[V]) != 0 ? times[V] / Sigma_a(V) : Real(0);
      for (ulong k = 0; k < size_3(Sigma_acc); ++k) {
        auto const index_k = order2 * (omega_end - omega_begin) * k;
        for (ulong j = 0; j < omega_end - omega_begin; ++j) {
          auto const index_j_k = index_k + j * order2;
          std::array<Complex, poly_n + 1> acc;
          acc.fill(Real(0));
          for (ulong m = 0; m < order; ++m) {
            auto const index_m_j_k = m * order + index_j_k;
            auto const index_m_V = m * order + index_V;
#pragma GCC unroll 4
#pragma GCC ivdep
            for (ulong l = 0; l < order; ++l) {
              auto const exp_K_c = exp_K[l + index_m_j_k];
              auto const Sigma_index = l + index_m_V;
              for (ulong p = 0; p <= p_max_S[order]; ++p) {
                acc[p] += exp_K_c * Sigma[Sigma_index][p];
              }  // for p in range(0, n_poly+1)
            }    // for l in range(0, bm_0_S)
          }      // for j in range(0, bm_0_S)
          for (ulong p = 0; p <= p_max_S[order]; ++p) {
            Sigma_acc(0, p, j, k, order) &= real(acc[p]) * one_Sigma_a_V;
            Sigma_acc(1, p, j, k, order) &= imag(acc[p]) * one_Sigma_a_V;
          }
        }  // for j in range(0, omega_end-omega_begin)
      }    // for k in range(0, size_y(acc))
    }
    if ((t_MC % accumulate_every) != 0)
      return;
    for (ulong S = 1; S < (1ul << order) - 1; ++S) {
      auto const bm_0_S = bitmap(0, S);
      auto const index_S = S * order2;
      auto const one_Sigma_a_S =
          std::abs(Sigma_a[S]) != 0 ? times[S] / Sigma_a(S) : Real(0);
      for (ulong k = 0; k < size_3(Sigma_acc); ++k) {
        auto const index_k = order2 * (omega_end - omega_begin) * k;
        for (ulong j = 0; j < omega_end - omega_begin; ++j) {
          auto const index_j_k = index_k + j * order2;
          std::array<Complex, poly_n + 1> acc;
          acc.fill(Real(0));
          for (ulong m = 0; m < bm_0_S; ++m) {
            auto const bm_m = bitmap(m + 1, S);
            auto const index_m_j_k = bm_m * order + index_j_k;
            auto const index_m_S = bm_m * order + index_S;
#pragma GCC unroll 4
#pragma GCC ivdep
            for (ulong l = 1; l <= bm_0_S; ++l) {
              auto const bm_l = bitmap(l, S);
              auto const exp_K_c = exp_K[bm_l + index_m_j_k];
              auto const Sigma_index = bm_l + index_m_S;
              for (ulong p = 0; p <= p_max_S[bm_0_S]; ++p) {
                acc[p] += exp_K_c * Sigma[Sigma_index][p];
              }  // for p in range(0, n_poly+1)
            }    // for l in range(0, bm_0_S)
          }      // for j in range(0, bm_0_S)
          for (ulong p = 0; p <= p_max_S[bm_0_S]; ++p) {
            Sigma_acc(0, p, j, k, bm_0_S) &= real(acc[p]) * one_Sigma_a_S;
            Sigma_acc(1, p, j, k, bm_0_S) &= imag(acc[p]) * one_Sigma_a_S;
          }  // for p in range(0, n_poly+1)
        }    // for j in range(0, omega_end-omega_begin)
      }      // for k in range(0, size_y(acc))
    }        // for S in range(0, (1ul<<order))
    for (ulong j = 0; j < size(Sigma_acc); ++j) {
      Sigma_acc[j].exec_deferred_add();
    }  // for j in sizerange(Sigma_acc)
  };

  std::string log_file_name{"../" + process_ID + "/log.log"};
  auto print_timers_r = [t_counter = (1. + Proba()) * print_interval, &timers,
                         &chosen_sets, log_file_name](ulong t_MC = 0) mutable {
    if (ffd::get<'T'>(timers)() > t_counter) {
      std::stringstream ss;
      auto print_timer = [&ss, &timers, t_MC](auto t,
                                              std::string const& x) mutable {
        ss << "t(" << x << ")= " << 1000 * 1000 * t() / t_MC << " \u03BCs   "
           << 100 * t() / ffd::get<'T'>(timers)() << " %\n";
        t.reset();
      };
      print_timer(ffd::get<'M'>(timers), "matrix  ");
      print_timer(ffd::get<'X'>(timers), "Xi      ");
      print_timer(ffd::get<'r'>(timers), "rho     ");
      print_timer(ffd::get<'S'>(timers), "Sigma   ");
      print_timer(ffd::get<'A'>(timers), "abs     ");
      print_timer(ffd::get<'H'>(timers), "heatbath");
      print_timer(ffd::get<'E'>(timers), "exp     ");
      print_timer(ffd::get<'K'>(timers), "acc_K   ");
      print_timer(ffd::get<'T'>(timers), "total   ");
      for (ulong j = 0; j < order + 1; ++j)
        ss << "card=" << j << ", t=" << chosen_sets[j] * 100. / t_MC << "%\n";
      if (t_MC != 0)
        ss << "after " << t_MC << " iterations "
           << "\n";
      if (verbose_stdout >= 2)
        std::cout << ss.str();
      if (verbose_file >= 2) {
        std::ofstream out_log(log_file_name);
        out_log << ss.str();
      }
      t_counter += (1. + Proba()) * print_interval;
    }
  };

  auto const one_factorial = [] {
    std::array<Real, order + 1> one_factorial;
    for (ulong u = 0; u < order + 1; ++u)
      one_factorial[u] = Real(1) / std::tgamma(u + 1);
    return one_factorial;
  }();
  auto print_accumulators_r = [p_max_n, p_max_S,
                               t_counter = (1. + Proba()) * print_interval,
                               &normalization_phase, &timers, &Sigma_acc,
                               &Un_acc, log_file_name, one_factorial,
                               Sigma_fnames = Sigma_fnames,
                               density_fname = density_fname]() mutable {
    if (ffd::get<'T'>(timers)() > t_counter && !normalization_phase) {
      auto const Sigma_norm = ENumber(Sigma_acc(0, 0, 0, 0, 0).mean_error());
      std::stringstream to_stdout, to_density_file, first_line_file;
      to_stdout << "Norm= " << std::setprecision(10) << to_string(Sigma_norm)
                << "\n";
      first_line_file << -1 << " " << std::setprecision(10) << -1 << " "
                      << Sigma_norm.value << " " << std::setprecision(3)
                      << Sigma_norm.error << "\n";
      to_density_file << 0 << " " << 0 << " " << n0 << " " << 0 << "\n";
      for (ulong u = 2; u < order + 1; ++u) {
        for (ulong p = (u == 2); p <= p_max_n[u]; ++p) {
          auto n_err = ENumber(Un_acc(p, u).mean_error());
          n_err *= 2 * one_factorial[u - 1];
          auto const n_err_norm = n_err / Sigma_norm;
          to_stdout << "n_{" << u - 1 << ", " << p
                    << "} = " << to_string(n_err_norm) << "\n";
          to_density_file << u - 1 << " " << p << " " << std::setprecision(8)
                          << n_err.value << " " << std::setprecision(3)
                          << n_err.error << "\n";
        }  // for p in range(0, poly_n+1)
      }    // for u in range(0, order+1)
      {
        std::ofstream density_f(density_fname, std::ofstream::trunc);
        density_f << first_line_file.str() << to_density_file.str();
      }
      for (ulong k = 0; k < size_3(Sigma_acc); ++k) {
        for (ulong j = 0; j < size_2(Sigma_acc); ++j) {
          std::stringstream to_file_r, to_file_i;
          for (ulong u = 2; u < order + 1; ++u) {
            for (ulong p = 0; p <= p_max_S[u]; ++p) {
              auto Sigma_r = ENumber(Sigma_acc(0, p, j, k, u).mean_error());
              auto Sigma_i = ENumber(Sigma_acc(1, p, j, k, u).mean_error());
              Sigma_r *= one_factorial[u];
              Sigma_i *= one_factorial[u];
              auto const Sigma_r_norm = Sigma_r / Sigma_norm;
              auto const Sigma_i_norm = Sigma_i / Sigma_norm;
              to_stdout << Sigma_fnames(0, j, k);
              to_stdout << ", " << u << " " << p << "= "
                        << to_string(Sigma_r_norm) << " +I*( "
                        << to_string(Sigma_i_norm) << " )\n";
              to_file_r << u << " " << p << " " << std::setprecision(8)
                        << Sigma_r.value << " " << std::setprecision(3)
                        << Sigma_r.error << "\n";
              to_file_i << u << " " << p << " " << std::setprecision(8)
                        << Sigma_i.value << " " << std::setprecision(3)
                        << Sigma_i.error << "\n";
            }  // for p in range(0, poly_n+1)
          }    // for u in range(1, order+1)
          {
            std::stringstream ss;
            ss << first_line_file.str() << to_file_r.str();
            std::ofstream self_f_r(Sigma_fnames(0, j, k), std::ofstream::trunc);
            self_f_r << ss.str();
          }
          {
            std::stringstream ss;
            ss << first_line_file.str() << to_file_i.str();
            std::ofstream self_f_i(Sigma_fnames(1, j, k), std::ofstream::trunc);
            self_f_i << ss.str();
          }
        }  // for j in range(0, size_x(Sigma_acc))
      }    // for k in range(0, size_y(Sigma_acc))
      if (verbose_stdout >= 3)
        std::cout << to_stdout.str();
      t_counter += (1. + Proba()) * print_interval;
    }
  };

  auto check_for_kill_signal = [t_counter =
                                    (1. + Proba()) * kill_signal_interval,
                                &timers, process_ID = process_ID]() mutable {
    if (ffd::get<'T'>(timers)() > t_counter) {
      if (fs::exists("kill") || fs::exists("../kill") ||
          fs::exists(process_ID + "/kill") || fs::exists("../../kill"))
        std::exit(0);
      t_counter += (1. + Proba()) * kill_signal_interval;
    }
  };

  for (ulong t_MC = 0;; ++t_MC) {
    ffd::get<'T'>(timers).ini();

    ffd::get<'M'>(timers).ini();
    fill_G0_matrix_r();
    ffd::get<'M'>(timers).fin();

    //    ffd::assert_special_case = true;
    ffd::get<'X'>(timers).ini();
    // ffd::sigmadet_hubbard::CDet_Xi_r(M, Sigma, bitmap);
    ffd::sigmadet_hubbard::CDet_Xi_poly_r(M0, Sigma, bitmap);
    ffd::get<'X'>(timers).fin();
    //    ffd::assert_special_case = false;
    ffd::get<'r'>(timers).ini();
    ffd::sigmadet_hubbard::CDet_rho_T_poly_r(M0, Sigma, rho_T, bitmap);
    ffd::get<'r'>(timers).fin();
    ffd::get<'S'>(timers).ini();
    ffd::sigmadet_hubbard::CDet_Sigma_r<order>(Sigma, rho_T, bitmap);
    ffd::get<'S'>(timers).fin();
#ifdef FFD_BUG_FLAG
    for (ulong S = 0; S < (1ul << order); ++S) {
      for (ulong j = 0; j < order; ++j) {
        for (ulong k = 0; k < order; ++k) {
          for (ulong p = 0; p < poly_n + 1; ++p) {
            std::cout << S << " (" << k << ", " << j << ") " << p << " "
                      << Sigma(k, j, S)[p] << "\n";
          }  // for p in range(0, poly_n+1)
        }    // for S in range(0, (1ul<<n))
      }      // for k in range(0, order)
    }        // for j in range(0, order)
    std::exit(0);
#endif

    ffd::get<'A'>(timers).ini();
    Sigma_a_r();
    ffd::get<'A'>(timers).fin();
    ffd::get<'H'>(timers).ini();
    auto const [times, X_new, chosen_set] =
        ffd::heat_bath_mc::HeatBathAll<order>(Sigma_a, proposer, X_i);
    ffd::get<'H'>(timers).fin();

    ++chosen_sets[ffd::core_math::popcount(chosen_set)];
    set_lambda_norm_r();

    ffd::get<'E'>(timers).ini();
    exp_K_r();
    ffd::get<'E'>(timers).fin();
    ffd::get<'K'>(timers).ini();
    Sigma_K_acc_r(times, t_MC);
    ffd::get<'K'>(timers).fin();

    X_i = X_new;

    ffd::get<'T'>(timers).fin();
    print_timers_r(t_MC + 1);
    print_accumulators_r();
    check_for_kill_signal();
  }

}  // main
