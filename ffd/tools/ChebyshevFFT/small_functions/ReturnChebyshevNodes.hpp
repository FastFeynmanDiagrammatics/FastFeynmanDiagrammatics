

namespace ffd::chebyshev_fft{

  std::vector<Real>

  ReturnChebyshevNodes(std::array<Real, 2> lower_upper_bound,
		       int log2order){
    unsigned long const num_nodes = 1<<log2order;
    std::vector<Real> nodes(num_nodes);
    using ffd::core_math::Pi;
    Real const pi_over_pow_2_order = Pi/num_nodes;


    auto [half_mean, half_diff] =
      return_half_mean_diff(lower_upper_bound);

    
    for(unsigned long j=0; j < num_nodes; ++j){
      nodes[j] = half_mean + half_diff*std::cos((j+0.5)*pi_over_pow_2_order);
    }

    
    return nodes;
  }

}//namespace
