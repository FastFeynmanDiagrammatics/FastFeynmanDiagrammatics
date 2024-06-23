namespace ffd::gauss_determinant{


  template<typename T>
  using WickMatrix = ffd::wick_matrix::WickMatrix<T>;

  

  template<typename Field>
  Field Determinant(WickMatrix<Field>& M_){
    using std::abs;
    if(size(M_) == 0)
      return 1;
    if(size(M_) == 1)
      return M_(0, 0);
    Field ret = 1;
    for(std::size_t j=0; j<size(M_); ++j){
      Real pivot_abs = abs(M_(j, j));
      std::size_t pivot_k = j;
      for(std::size_t k=j+1; k<size(M_); ++k){
	if(abs(M_(j, k)) > pivot_abs){
	  pivot_abs = abs(M_(j, k));
	  pivot_k = k;
	}
      }
      if( pivot_abs < 10000*std::numeric_limits<Real>::min() ){
	return Field(0.);
      }
      if(pivot_k != j){
	for(std::size_t j2=j; j2<size(M_); ++j2){
	  Field temp_mat = M_(j2, j);
	  M_(j2, j, 'a') = M_(j2, pivot_k);
	  M_(j2, pivot_k, 'a') = temp_mat;
	}
	ret = ret*(-1);
      }
      Field temp_ret = M_(j,j)*ret;
      ret = temp_ret;
      for(std::size_t j2=j+1; j2<size(M_); ++j2){
	Field temp_mat = M_(j2, j)/M_(j,j);
	M_(j2, j, 'a') = temp_mat;
      }
      for(std::size_t k=j+1; k < size(M_); ++k){
	for(std::size_t j2=j+1; j2 < size(M_); ++j2){
	  M_(j2, k, 'a') -= M_(j,k)*M_(j2,j);
	}
      }
    }
    return ret;
  }


}//namespace ffd::gauss_determinant

