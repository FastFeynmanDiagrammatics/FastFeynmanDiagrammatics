namespace ffd::ranked_zeta::unit_test{

	void abs(){

  	  using NilpotentPolynomial = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Full>;
  	  using EvenNilPoly = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Even>;
  	  using namespace user_space;

	  const int order = 10;
	  NilpotentPolynomial P1(order);
	  for(int j=0; j < 1<<order; ++j){
		  P1[j] = std::sin(0.7+4.4*j);
		}

	EvenNilPoly P1_e(order);
	for(int j=0; j < 1<<order; ++j){
		if (!__builtin_parity(j)){
			P1_e[j] = std::sin(0.7+4.4*j);
		}
	  }

	  RankedZeta Z1 = P1;
	  EvenZeta Z1_e = P1_e;

	  P1 = abs(P1);
  	  P1_e = abs(P1_e);

	  Z1=abs(Z1);
	  Z1_e=abs(Z1_e);

	  NilPoly P1_new=Z1;
  	  EvenNilPoly P1_e_new=Z1_e;

	  for(int j=0; j < 1<<order; ++j){
	   assert(std::abs(P1[j] - P1_new[j]) < 100*std::numeric_limits<long double>::epsilon() );
	   assert(std::abs(P1_e[j] - P1_e_new[j]) < 100*std::numeric_limits<long double>::epsilon() );
	  }
    }

}//namespace
