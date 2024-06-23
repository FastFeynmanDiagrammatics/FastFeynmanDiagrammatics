

namespace ffd::user_space::imaginary_time_lattice_propagator{

  template<int d,
	   int n,
	   typename T>
  void
  get_linear_sizes(T const& x,
		   std::vector<int>& linear_sizes){
    if constexpr(d > 0){
	auto const& low_upp = ffd::get<n>(x).LowerUpperBound.value();
	
	linear_sizes[n] = low_upp[1]-low_upp[0];
	if constexpr( n < d-1){
	    get_linear_sizes<d, n+1>(x, linear_sizes);
	  }
      }
  }


}//namespace
