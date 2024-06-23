namespace ffd::determinant{

  template<typename matrix_t>

  typename std::decay<decltype(std::declval<matrix_t>()[0])>::type
  
  Determinant(matrix_t&& M_){
    using std::abs;
    using element_t = typename std::decay<decltype(std::declval<matrix_t>()[0])>::type;

    
    std::size_t const lin_size = ffd::core_math::sqrt_int( size(M_) );
    assert(lin_size*lin_size == size(M_) );

    
    switch(lin_size){
    case 0ul: 
      return element_t(1.);
    case 1ul:
      return M_[0];
    case 2ul:
      return M_[0]*M_[3] - M_[1]*M_[2];
    case 3ul:
      return
	M_[0] * (M_[4]*M_[8] - M_[5]*M_[7]) -
	M_[1] * (M_[3]*M_[8] - M_[5]*M_[6]) +
	M_[2] * (M_[3]*M_[7] - M_[4]*M_[6]);
    case 4ul:
      return
	+ (M_[0]*M_[5] - M_[4]*M_[1]) * (M_[10]*M_[15] - M_[11]*M_[14])
	- (M_[0]*M_[6] - M_[4]*M_[2]) * (M_[9]*M_[15] - M_[11]*M_[13])
	+ (M_[0]*M_[7] - M_[4]*M_[3]) * (M_[9]*M_[14] - M_[10]*M_[13])
	+ (M_[1]*M_[6] - M_[5]*M_[2]) * (M_[8]*M_[15] - M_[11]*M_[12])
	- (M_[1]*M_[7] - M_[5]*M_[3]) * (M_[8]*M_[14] - M_[10]*M_[12])
	+ (M_[2]*M_[7] - M_[6]*M_[3]) * (M_[8]*M_[13] - M_[9]*M_[12]);
    default:
      element_t det = 1.;

      
      for(std::size_t j=0; j<lin_size-1; ++j){
	Real pivot_abs = abs( M_[j+lin_size*j] );
	std::size_t pivot_k = j;
	for(std::size_t k=j+1; k<lin_size; ++k){
	  if(  abs( M_[j+k*lin_size] )  >
	       pivot_abs  ){
	    pivot_abs = abs( M_[j+k*lin_size] );
	    pivot_k = k;
	  }
	}
	if( pivot_abs < 10000*std::numeric_limits<Real>::min() ){
	  return element_t(0.);
	}


	//SWAP ROWS
	if(pivot_k != j){
	  for(std::size_t j2=j; j2<lin_size; ++j2){
	    element_t temp_mat = M_[j2 + j*lin_size];
	    M_[j2+j*lin_size] = M_[j2+pivot_k*lin_size];
	    M_[j2+pivot_k*lin_size] = temp_mat;
	  }
	  det = -det;
	}



	//GAUSSIAN ELIMINATION
	element_t temp_det = M_[j+j*lin_size]*det;
	det = temp_det;
	element_t const one_over_pivot = element_t(1.)/M_[j+j*lin_size];
	for(std::size_t j2=j+1; j2<lin_size; ++j2){
	  M_[j2+ j*lin_size] = M_[j2+j*lin_size]*one_over_pivot;
	}
	for(std::size_t k=j+1; k < lin_size; ++k){
	  for(std::size_t j2=j+1; j2 < lin_size; ++j2){
	    M_[j2+ k*lin_size] -= M_[j+k*lin_size]*M_[j2+j*lin_size];
	  }
	}
      }

    
      return det*M_[lin_size*lin_size-1];
    }
  }
  


  }//namespace
