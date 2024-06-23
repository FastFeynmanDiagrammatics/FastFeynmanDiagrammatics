

namespace ffd::user_space::imaginary_time_lattice_propagator::unit_test{

  template<int d, bool is_fermion = true>
  struct HypercubicPropagatorExact{
    static constexpr auto I =  ffd::user_space::I;
    
    
    static constexpr Real zeta = 1 - 2*(is_fermion);

    
    std::array<int, d> L;
    std::array<Real, 2> Mu;
    std::array<std::array<Real, d>, 2> t_spin_j;
    Real beta;

    HypercubicPropagatorExact() = default;

    
    long n_omega_max = 1 << 12;

    
    Real xi_k(std::array<Real, d> k,
	      int spin) const{
      Real ret = -Mu[(1+spin)/2];
      for( int j: ffd::vector_range::Range(d) ){
	ret += -2*t_spin_j[(1+spin)/2][j]*cos(k[j]);
      }
      return ret;
    }
    
    
    
     

    Real
    operator()(Real tau,
	       std::array<int, d> r,
	       int spin) const{
      using ffd::user_space::expI;
      using ffd::vector_range::Range;
      using ffd::core_math::Pi;
      
      
      Complex propagator = 0.;

      
      std::array<std::vector<Real>, d> k_ranges;
      for( int j: Range(d) ){
	for( int n: Range( L[j] ) ){
	  k_ranges[j].push_back( 2*n*ffd::core_math::Pi / L[j] );
	}
      }
      
      
      for( auto k: ffd::cartesian_product_array::
	     CartesianProductArray<Real, d>(k_ranges) ){
	for( int n_omega: Range(-n_omega_max, n_omega_max) ){


	  Real k_dot_r = 0;
	  for( int j: Range(d) ){
	    k_dot_r += k[j]*r[j];
	  }
	  
	  
	  Real omega_n = 2*Pi*(n_omega+.5)/beta;


	  propagator += expI( k_dot_r - omega_n*tau )/
	    (I*omega_n - xi_k(k, spin));
	  
	}
      }


      Real ret = std::real(propagator)/beta;
      for( int j: Range(d) ){
	ret /= L[j];
      }
      
      
      return ret;

    }

    
    Real
    operator()(bool,
	       Real tau,
	       std::array<int, d> r,
	       int spin) const{
      using ffd::user_space::expI;
      using ffd::vector_range::Range;
      using ffd::core_math::Pi;
      
      
      Complex propagator = 0.;

      
      std::array<std::vector<Real>, d> k_ranges;
      for( int j: Range(d) ){
	for( int n: Range( L[j] ) ){
	  k_ranges[j].push_back( 2*n*ffd::core_math::Pi / L[j] );
	}
      }
      
      
      for( auto k: ffd::cartesian_product_array::
	     CartesianProductArray<Real, d>(k_ranges) ){
	

	Real k_dot_r = 0;
	  for( int j: Range(d) ){
	    k_dot_r += k[j]*r[j];
	  }


	  Real E = xi_k(k, spin);
	  

	  propagator += -expI( k_dot_r )*exp( - tau*E )/
	    (1-zeta*exp(-beta*E));
	  
      }

      
      Real ret = std::real(propagator);
      for( int j: Range(d) ){
	ret /= L[j];
      }
      return ret;
    }
    
    
  };

}//namespace
