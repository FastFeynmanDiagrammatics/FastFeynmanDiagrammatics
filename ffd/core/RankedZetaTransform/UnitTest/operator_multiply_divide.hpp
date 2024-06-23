namespace ffd::ranked_zeta::unit_test{

	void operator_multiply_divide(){
	  using NilpotentPolynomial = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Full>;
	  using EvenNilPoly = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Even>;
	  using namespace user_space;
    using std::size;

	  const int order = 10;
	  NilpotentPolynomial P1(order), P2(order);
	  for(int j=0; j < 1<<order; ++j){
		  P1[j] = std::sin(0.7+4.4*j);
		  P2[j] = std::sin(1.5+2.3*j);
		}

  	auto P_mul = P1*P2;
  	auto P_div = P1/P2;

  	RankedZeta Z1 = P1;
  	RankedZeta Z2 = P2;

  	auto Z_mul = Z1*Z2;
  	auto Z_div = Z1/Z2;

    auto Z1_ = Z1;
    Z1*=Z2;
    Z1_/=Z2;

  	NilpotentPolynomial P_mul_new = Z_mul;
  	NilpotentPolynomial P_div_new = Z_div;

    NilpotentPolynomial P_mul_alt = Z1;
    NilpotentPolynomial P_div_alt = Z1_;


  	for(int j=0; j < 1<<order; ++j){
      // std::cerr << j << " " <<  P_mul[j] << " " << P_mul_new[j] << " " << std::abs(P_mul[j] - P_mul_new[j]) <<  std::endl;
      // std::cerr << j << " " <<  P_div[j] << " " << P_div_new[j] << " " << std::abs(P_div[j] - P_div_new[j]) <<  std::endl;
  	  assert(std::abs((P_mul[j] - P_mul_new[j])/P_mul[j])  < 1.e-10);
      assert(std::abs((P_div[j] - P_div_new[j])/P_div[j])  < 1.e-10);
      assert(std::abs((P_mul[j] - P_mul_alt[j])/P_mul[j])  < 1.e-10);
      assert(std::abs((P_div[j] - P_div_alt[j])/P_div[j])  < 1.e-10);
  	}

    EvenNilPoly P1_e(order), P2_e(order);
  	for(int j=0; j < 1<<order; ++j){
  		if (!__builtin_parity(j)){
  			P1_e[j] = std::sin(0.7+4.4*j);
  			P2_e[j] = std::sin(1.5+2.3*j);
  		}
  	}

  	auto P_mul_e = P1_e*P2_e;
  	auto P_div_e = P1_e/P2_e;

  	RankedZeta Z1_e = P1_e;
  	RankedZeta Z2_e = P2_e;

  	auto Z_mul_e = Z1_e*Z2_e;
  	auto Z_div_e = Z1_e/Z2_e;

    auto Z1_e_ = Z1_e;
    Z1_e*=Z2_e;
    Z1_e_/=Z2_e;

  	EvenNilPoly P_mul_new_e = Z_mul_e;
  	EvenNilPoly P_div_new_e = Z_div_e;

    EvenNilPoly P_mul_alt_e = Z1_e;
    EvenNilPoly P_div_alt_e = Z1_e_;

  	for(int j=0; j < 1<<order; ++j){
      // std::cerr << j << " " <<  P_mul_e[j] << " " << P_mul_new_e[j] << " " << std::abs(P_mul_e[j] - P_mul_new_e[j]) <<  std::endl;
      // std::cerr << j << " " <<  P_div_e[j] << " " << P_div_new_e[j] << " " << std::abs(P_div_e[j] - P_div_new_e[j]) <<  std::endl;
      if (P_mul_e[j]!=0. || P_mul_new_e[j]!=0.) assert(std::abs((P_mul_e[j] - P_mul_new_e[j])/P_mul_e[j])  < 1.e-10);
      if (P_div_e[j]!=0. || P_div_new_e[j]!=0.) assert(std::abs((P_div_e[j] - P_div_new_e[j])/P_div_e[j])  < 1.e-10);
      if (P_mul_e[j]!=0. || P_mul_new_e[j]!=0.) assert(std::abs((P_mul_e[j] - P_mul_alt_e[j])/P_mul_e[j])  < 1.e-10);
      if (P_div_e[j]!=0. || P_div_new_e[j]!=0.) assert(std::abs((P_div_e[j] - P_div_alt_e[j])/P_div_e[j])  < 1.e-10);
    }

    // const Real x = 2.;
    const Real x = 1.23456789123123131231231242l;
    // FullZeta Zi = x;
    // FullZeta Zj = x;
    // FullZeta Zk = x;
    // assert(size(Zi) == 1);

    auto Pii = P2*x;
    auto Pjj = x/P2;
    auto Pkk = P2/x;

    auto Pii_e = P2_e*x;
    auto Pjj_e = x/P2_e;
    auto Pkk_e = P2_e/x;

    auto Zi = x*Z2;
    auto Zj = x/Z2;
    auto Zk = Z2/x;

    auto Zi_e = x*Z2_e;
    auto Zj_e = x/Z2_e;
    auto Zk_e = Z2_e/x;

    NilpotentPolynomial Pi = Zi;
    NilpotentPolynomial Pj = Zj;
    NilpotentPolynomial Pk = Zk;

    EvenNilPoly Pi_e = Zi_e;
    EvenNilPoly Pj_e = Zj_e;
    EvenNilPoly Pk_e = Zk_e;

    for(BinaryInt j=0; j < size(Pi); ++j){
      assert(std::abs(Pi[j] - Pii[j]) < std::abs(Pii[j])*1.e-10);
      assert(std::abs(Pj[j] - Pjj[j]) < std::abs(Pjj[j])*1.e-10);
      assert(std::abs(Pk[j] - Pkk[j]) < std::abs(Pkk[j])*1.e-10);
      if (Pi_e[j]!=0. || Pii_e[j]!=0.) assert(std::abs(Pi_e[j] - Pii_e[j]) < std::abs(Pii_e[j])*1.e-10);
      if (Pj_e[j]!=0. || Pjj_e[j]!=0.) assert(std::abs(Pj_e[j] - Pjj_e[j]) < std::abs(Pjj_e[j])*1.e-10);
      if (Pk_e[j]!=0. || Pkk_e[j]!=0.) assert(std::abs(Pk_e[j] - Pkk_e[j]) < std::abs(Pkk_e[j])*1.e-10);
    }

  }

}//namespace
