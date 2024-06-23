namespace ffd::normal_principal_blocks{

  template<typename T>
  void
  pad_matrix(std::vector<T>& M,
	     int p = -1){
    std::size_t const l_size = ffd::core_math::sqrt_int(size(M));
    if(p == -1) p = l_size;
    std::vector<T> M_p((l_size+1)*(l_size+1));
    for(uint j=0; j<l_size+1; ++j){
      uint j_M = j < p ? j : j - 1;
      for(uint k=0; k<l_size+1; ++k){
	if(j != p && k != p){
	  uint k_M = k < p ? k : k - 1;
	  M_p[k+(l_size+1)*j] = M[k_M+l_size*j_M];
	}else{
	  M_p[k+(l_size+1)*j] = Real((j == k));
	}
      }
    }
    M = M_p;
  }


}//namespace
