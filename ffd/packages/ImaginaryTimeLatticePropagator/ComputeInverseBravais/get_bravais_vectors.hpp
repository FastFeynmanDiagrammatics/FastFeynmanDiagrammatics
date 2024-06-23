

namespace ffd::user_space::imaginary_time_lattice_propagator{

  template<int d,
	   int n,
	   typename T>
  void
  get_bravais_vectors(T x,
		      std::vector<Real>& bravais_vectors){
    if constexpr( d != 0){
	std::array<Real, d> bravais = ffd::get<n>(x).RealSpaceVectors[0];
	for(int j=0; j < d; ++j){
	  bravais_vectors[j*d + n] = bravais[j];
	}
	if constexpr( n < d-1){
	    get_bravais_vectors<d, n+1>(x, bravais_vectors);
	  }
      }
  }




}//namespace
