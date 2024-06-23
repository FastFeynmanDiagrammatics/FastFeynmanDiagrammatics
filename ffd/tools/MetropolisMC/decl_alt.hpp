namespace ffd::metropolis_mc {

// we sample with w = |f(x)| + lambda*g(x)
// we store the observables
// <O_1> = <f(x)/[|f(x)| + lambda*g(x)]>
// and
// <O_2> = <|f(x)|/[|f(x)| + lambda*g(x)]>
// where
// \int g(x) = 1
// where g is the product of conditional probabilities of creating X[i] from the origin.
// We have
// int f(x) / int |f(x)| = <O_1>/<O_2>
// and
// <O_2> = int |f(x)|/(int |f(x)|+lambda*1)
// which yields
// int |f(x)| = (lambda*<O_2>)/(1-<O_2>)
// so that
// int f(x) = <O_1>/<O_2> * (lambda*<O_2>)/(1-<O_2>).

  template <class Field_f = Real,
	          class Field_g = Real,
	          class Field_l = Real>
  auto MC_step_hybrid(
		      std::array<Field_f, 2>& f,
		      std::array<Field_g, 2>& g,
		      Field_l const lambda
		      ) {
    using std::abs;

    // compute weights
    std::array<std::decay_t<decltype(abs(g[0]))>, 2> w;
    for (auto s: {0, 1} ) {
      w[s] = abs(f[s]) + lambda*abs(g[s]);
    }

    // perform metropolis update with some probability
    const auto ratio = w[1] / w[0];
    auto c = (ratio > random_distributions::Proba());

    f[0] = f[c];
    g[0] = g[c];

    return std::make_pair(c, w[c]);
  }

} // namespace ffd::metropolis_mc
