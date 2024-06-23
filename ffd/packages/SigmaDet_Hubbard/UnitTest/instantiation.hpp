namespace ffd::sigmadet_hubbard::unit_test{

  template<int n>
  void instantiation(){
    using nilpoly_t = ffd::nilpotent_polynomial::NilpotentPolynomial<Real>;
    BinaryInt constexpr two_n = (1u<<n);
    std::array<nilpoly_t, n*n> G0;
    


    for(int j=0; j< n; ++j){
      for(int k = 0; k < n; ++k){
	G0[j*n+k] = nilpoly_t(n);
	for( BinaryInt S = 0; S < two_n; ++S){
	  if( j != k){
	    G0[j*n+k][S] = ffd::core_math::Factorial(ffd::set_theory::CardinalitySet(S))*(1-2*(j>k));
	  }
	}
	// std::cerr<<j<<": "<<G0[j]<<std::endl;
      }
    }


    [[maybe_unused]] auto Xi = compute_Xi(G0);
    // for(int j=0; j< n; ++j){
    //   for(int k=0; k < n; ++k){
    // 	std::cerr<<"Xi["<<j<<", "<<k<<"] = "<<Xi[j*n+k]<<std::endl;
    //   }
    // }

    
    [[maybe_unused]] auto rho = compute_rho(G0, Xi);
    // for(int j=0; j< n; ++j){
    //   for(int k=0; k < n; ++k){
    // 	std::cerr<<"rho["<<j<<", "<<k<<"] = "<<rho[j*n+k]<<std::endl;
    //   }
    // }


    [[maybe_unused]] auto Sigma = compute_Sigma(Xi, rho);
    // for(int j=0; j< n; ++j){
    //   for(int k=0; k < n; ++k){
    // 	std::cerr<<"Sigma["<<j<<", "<<k<<"] = "<<Sigma[j*n+k]<<std::endl;
    //   }
    // }

  }

}//namespace
