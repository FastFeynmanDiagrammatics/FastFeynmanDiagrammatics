namespace ffd::user_space{

  template<int Order,
	   typename accumulators_t,
	   typename array_times_t>
    
  void
  add_measurements_to_accumulators(accumulators_t& accumulators,
				   array_times_t const& times){
    using std::abs, std::sqrt;
    BinaryInt const two_Order = (1<<Order);

    
    std::array<Real, Order+1> one_over_factorial_binomial;
    for(std::size_t j=0; j<=Order; ++j){
      one_over_factorial_binomial[j] = 1./
	(
	 ffd::core_math::Factorial(j)*
	 ffd::core_math::BinomialCoefficient(Order, j)
	 );
    }

    
    std::array<Real, Order+1> acc;
    acc.fill(0.);

    
    for(BinaryInt S = 0; S < two_Order; ++S){
      acc[__builtin_popcount(S)] += times[S];
    }
    
    for(int j=0; j<=Order; ++j){
      accumulators[j] += acc[j]*one_over_factorial_binomial[j];
    }
    
  }

}//namespace
