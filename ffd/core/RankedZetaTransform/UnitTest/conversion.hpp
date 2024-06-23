namespace ffd::ranked_zeta::unit_test{

  void conversion(){

    using NilpotentPolynomial = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Full>;
	using EvenNilPoly = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Even>;
	using namespace user_space;

	const int order = 10;
    NilpotentPolynomial P(order);
	for(int j=0; j < 1<<order; ++j){
      P[j] = std::sin(1+j);
    }

    RankedZeta Z = P;

    assert(Z.NotZero);
	NilpotentPolynomial P_new = Z;

	for(int j=0; j < 1<<order; ++j){
       assert(std::abs(P[j] - P_new[j]) < 100*std::numeric_limits<long double>::epsilon() );
    }

	EvenNilPoly P_even(order);
	for(int j=0; j < 1<<order; ++j){
	  if (!__builtin_parity(j)){
	  	P_even[j] = std::sin(1+j);
	  }
	}

	EvenZeta Z_even = P_even;
	assert(Z.NotZero);

	// for(int j=0; j < 1<<order; ++j){
	// 	std::cerr<<j << "  ";
	// 	for (int k=0; k<Z_even.rank_size+1; ++k){
	// 		std::cerr << Z_even[j+k*Z_even.cardinality]<< " ";
	// 	}
	// 	std::cerr << std::endl;
	// }

	EvenNilPoly P_even_new = Z_even;

	for(int j=0; j < 1<<order; ++j){
	   assert(std::abs(P_even[j] - P_even_new[j]) < 100*std::numeric_limits<long double>::epsilon() );
	}

  }

}//namespace
