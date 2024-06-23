namespace ffd::metropolis_mc::unit_test {

void metropolisMC_test() {
  //  using namespace ffd::user_space;

  double norm_tot = 0.;
  double phys_tot = 0.;
  double g_tot = 0;
  const double O = 0.;
  double V = O;
  auto V_old = std::array<double, 1>{V};
  const double sigma = 0.5;
  const double mu = O;
  auto lambda = 15.0; // initial guess
  const auto N_steps = (1 << 23);
  const auto Therm_steps = (1 << 20);
  const auto Therm_rounds = 10;
  const auto target_factor = 10.;

  // f(x) = x*\theta(-x+2)*\theta(x) + 0.5*x*\theta(x+2)*\theta(-x)
  // \int f(x) = 1
  // \int |f(x)| = 3
  auto calc_w = [](auto x) {
    auto x_lim = 2.;
    if (x[0] < x_lim && x[0] > 0.) {
      return 1.0*x[0];
    } else if (x[0] < 0. && x[0] > -x_lim) {
      return 0.5*x[0];
    }
    return 0.;
  };

  // proposer generates a new position from normal distribution
  auto prop = [sigma, mu](auto x) {
    std::normal_distribution<Real> ND(mu, sigma);
    return x + ND(random_distributions::RNGen);
  };

  // computes conditional probability of generating x from normal distribution around mu
  auto cond_prop = [sigma](auto x, auto mu) {
    using core_math::Pi;
    auto term = (x - mu) / sigma;
    return exp(-0.5 * term * term) / (sigma * sqrt(2 * Pi));
  };

  // initialise
  auto f_old = calc_w(V_old);
  auto g_old = cond_prop(V_old[0], O);
  auto w_old = abs(f_old) + lambda * g_old;
  auto o_1_old = f_old      /(abs(f_old) + lambda*g_old);
  auto o_2_old = abs(f_old) /(abs(f_old) + lambda*g_old);

  // Thermalisation for lambda
  for (uint round = 0; round < Therm_rounds; ++round){
    for (uint steps = 0; steps < Therm_steps; ++steps) {
      // Metropolis Monte Carlo step
      auto [V_new, w_new, o_1_new, o_2_new, g_new] =
          MC_step<1, double>(V_old, w_old, o_1_old, o_2_old, g_old, prop, cond_prop, calc_w, O, lambda);
      // update
      V_old = V_new;
      w_old = w_new;
      o_1_old = o_1_new;
      o_2_old = o_2_new;
      g_old = g_new;
      //accumulate
      norm_tot += o_2_old;
      g_tot += g_old;
    }
    g_tot /= (double)Therm_steps;
    norm_tot /= (double)Therm_steps;
    auto norm_ff = (lambda * norm_tot)/(1. - norm_tot);
    auto lambda_new = norm_ff/(target_factor*g_tot);

    std::cout << "after integration lambda*g(x)   = " << g_tot*lambda << ", should be: " << norm_ff/target_factor << std::endl;
    std::cout << "new lambda is: " << lambda_new << std::endl;

    lambda = lambda_new;

    norm_tot = 0.;
    g_tot = 0.;
  }

  // do MC
  for (uint steps = 0; steps < N_steps; ++steps) {
    // Metropolis Monte Carlo step
    auto [V_new, w_new, o_1_new, o_2_new, g_new] =
        MC_step<1, double>(V_old, w_old, o_1_old, o_2_old, g_old, prop, cond_prop, calc_w, O, lambda);
    // update
    V_old = V_new;
    w_old = w_new;
    o_1_old = o_1_new;
    o_2_old = o_2_new;
    g_old = g_new;

    //accumulate
    phys_tot += o_1_old;
    norm_tot += o_2_old;
    g_tot += g_old;
  }

  phys_tot /= (double)N_steps;
  norm_tot /= (double)N_steps;
  g_tot /= (double)N_steps;
  auto const ratio = phys_tot / norm_tot;
  auto const normf = (lambda * norm_tot)/(1. - norm_tot);

  std::cout << "MC finished after " << N_steps << " steps." << std::endl;
  std::cout << "after integration g(x)          = " << g_tot << ", should be: " << 1.0 / (3.0 + lambda)  << std::endl;
  std::cout << "after integration f(x)/|f(x)|   = " << ratio << ", should be: 0.3333" << std::endl;
  std::cout << "after integration |f(x)|        = " << normf << ", should be: 3.0000" << std::endl;
  std::cout << "after integration f(x)          = " << ratio*normf << ", should be: 1.0000" << std::endl;
  std::cout << "after integration lambda*g(x)   = " << g_tot*lambda << ", should be: " << normf/target_factor << std::endl;

  auto sum = 0.;
  for (int i=-100; i<=100; ++i){
    sum += cond_prop(double(i),0.);
  }
  std::cout << sum << std::endl;

}
} // namespace ffd::metropolis_mc::unit_test
