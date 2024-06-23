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
template <int order, typename Field = double, typename seed_t = int,
          typename cond_t = int, typename coord_t = int,
          typename coord_O_t = int, typename weight_t = int>
auto MC_step(coord_t const &X_old, 
	     Field const w_old, 
	     Field const o_1_old, 
	     Field const o_2_old, 
	     Field const g_old, //remove later
	     seed_t prop,
             cond_t cond_prop, 
	     weight_t const &calc_w, 
	     coord_O_t const &O,
             Field const lambda) {

  // generate one new vertex from one randomly chosen old one	
  auto X_new = X_old;
  const int ind = ffd::random_distributions::RandomInRange(0, order);
  X_new[ind] = prop(X_old[ind]);
  
  // compute new weight (could be made faster)
  const auto f = calc_w(X_new);
  //auto f = 0.;
  auto g = 1.;
  for (uint j = 0; j < order; ++j) {
    g *= cond_prop(X_new[j], O);
  }
 
  const auto w_new = abs(f) + lambda * g;

  // perform metropolis update with some probability
  const auto ratio = w_new / w_old;
  const auto rn = random_distributions::Proba();
  if (ratio > rn) {  //accept
    const auto o_1_new =     f  / w_new;
    const auto o_2_new = abs(f) / w_new;    
    const auto g_new = g / w_new;
    return std::make_tuple(X_new, w_new, o_1_new, o_2_new, g_new);
  }
  // reject
  return std::make_tuple(X_old, w_old, o_1_old, o_2_old, g_old);
}

} // namespace ffd::metropolis_mc
