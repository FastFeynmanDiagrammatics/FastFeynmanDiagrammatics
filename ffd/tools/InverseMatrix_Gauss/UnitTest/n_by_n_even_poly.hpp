namespace ffd::inverse_matrix_gauss::unit_test{

  template<int n, int m>
  
  auto
  n_by_n_even_poly(Real factor_eps_real = 1ul<<25){
    using namespace ffd::user_space;
    bool IsOk = true;
    Real maximal_error = 0;
    using std::abs;
    using nilpoly_t = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Even>;
    long double const eps_real = factor_eps_real*std::numeric_limits<Real>::epsilon();
    std::array<nilpoly_t, n*n> M;
    for( int j=0 ; j< n*n; ++j){
      M[j] = nilpoly_t(m);
      for( BinaryInt S=0; S < (1<<m); ++S){
	M[j][S] = (2*ffd::user_space::Proba()-1.)*
	  ffd::core_math::Factorial(ffd::set_theory::CardinalitySet(S));
      }
    }
    auto M_copy = M;
    
    
    auto [inverse, det] = Inverse_and_Determinant(M);
    decltype(M) identity;
    for( int j=0; j < n; ++j ){
      for( int k=0; k < n; ++k ){
	for( int l=0; l < n; ++l){
	  identity[j*n+k] += inverse[j*n+l]*M_copy[l*n+k];
	}
	// std::cerr<<j<<" "<<k<<": "<<std::setprecision(1)<<identity[j*n+k]-(j==k)*1.<<"\n";
	
	
	for( BinaryInt S = 0; S < (1<<m); ++S){
	  if(__builtin_parity(S)==0){
	    Real const abs_error = abs(identity[j*n+k][S] - (S==0)*(j==k)*1.);
	    if(abs_error > maximal_error){
	      maximal_error = abs_error;
	    }
	    IsOk = IsOk &&
	      abs_error
	      < eps_real;
	  }
	}
      }
    }
    

    auto det2 = ffd::determinant::Determinant(M_copy);
    // std::cerr<<std::setprecision(20)<<"det = "<<det2<<" "<<det<<std::endl;
    Real maximal_error_det = 0;
    for( BinaryInt S = 0; S < (1<<m); ++S){
      Real const abs_error = abs(det[S]-det2[S]);
      if(abs_error > maximal_error_det){
	maximal_error_det = abs_error;
      }
      IsOk = IsOk && abs_error < eps_real;
    }


    return std::make_tuple(IsOk, maximal_error, maximal_error_det);
  }

}//namespace
