

namespace ffd::user_space::imaginary_time_lattice_propagator{

  std::vector<Real>
  ReturnWaveVectorRange(int L, std::optional<Real> alpha_twisted_boundary_conditions){
    using ffd::vector_range::Range;
    using ffd::core_math::Pi;

    std::vector<Real> wave_vector_range(L, 0.);

    int iter = 0;
    for( int u: Range(-((L-1)/2), L/2 + 1) ){
      
      wave_vector_range[iter] = 2*u*Pi/L;
      
      if( alpha_twisted_boundary_conditions.has_value() ){
	wave_vector_range[iter] += alpha_twisted_boundary_conditions.value()/L;
      }
      
      ++iter;
      
    }
    assert(( iter == L ));

    return wave_vector_range;
  }

}//namespace
