namespace ffd::gauss_pfaffian{

  template<typename antisymmetric_matrix_t>

  typename std::decay<decltype(std::declval<antisymmetric_matrix_t>()[0])>::type
  Pfaffian(antisymmetric_matrix_t&& A){
    using std::abs;

    
    int const n = return_n( size(A) );
    if( n == 0 ){ return 1.; }

    
    typename std::decay<decltype(A[0])>::type ret = 1.;

    
    for( int r = 0; r < n - 1; ++r){
      int p = 2*r + 1;
      Real pivot_max = abs(  A[ lin(2*r, p, n) ]  );
      for( int c = 2*r+2; c < 2*n; ++c){
	Real c_val = abs(  A[ lin(2*r, c, n) ]  ); 
	if( c_val > pivot_max ){
	  p = c;
	  pivot_max = c_val; 
	}
      }


      if( pivot_max <
	  100000*std::numeric_limits<Real>::min() ){
	return 0.;
      }


      if( p != 2*r + 1 ){
	std::swap( A[ lin(2*r, 2*r+1, n) ],
		   A[ lin(2*r, p, n) ] );


	A[ lin(2*r+1, p, n) ] *= -1.;


	for( int j = 2*r + 2; j < p; ++j){
	  std::swap( A[ lin(2*r+1, j, n) ],
		     A[ lin(j, p, n) ] );

	  
	  A[ lin(2*r+1, j, n) ] *= -1.;
	  A[ lin(j, p, n) ] *= -1.;
	}


	for( int j = p+1; j < 2*n; ++j){
	  std::swap( A[ lin(2*r+1, j, n) ],
		     A[ lin(p, j, n) ] );
	}


	ret *= -1.;
      }


      ret *= A[ lin(2*r, 2*r+1, n) ];


      for( int c = 2*r + 2; c < 2*n; ++c){
	auto lambda = A[ lin(2*r, c, n) ] / A[ lin(2*r, 2*r+1, n) ];

	
	for( int j = 2*r+2; j < c; ++j){
	  A[ lin(j, c, n) ] += lambda*A[ lin(2*r+1, j, n) ]; 
	}


	for( int j = c+1; j < 2*n; ++j){
	  A[ lin(c, j, n) ] -= lambda*A[ lin(2*r+1, j, n) ]; 
	}

	
      }
    }


    return ret*A[ lin(2*n-2, 2*n-1, n) ];
  }

}//namespace
