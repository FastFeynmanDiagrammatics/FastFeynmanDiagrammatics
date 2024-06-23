

namespace ffd::nilpotent_polynomial::unit_test{

  template<SetT SetType=Full>
  void test_divide(int order, int factor, int omega1, int omega2, int delta1, int delta2){
    NilpotentPolynomial<long double, SetType> P(order), Q(order), AP(order), AQ(order);
    for(int j=0; j < 1<<order; ++j){
	  if (SetType==Full || (SetType==Even && __builtin_parity(j)==0)){	// || (SetType==Odd && __builtin_parity(j)==1)
        P[j] = sin(omega1*j+delta1);
        AP[j] = std::abs(P[j]);
        Q[j] = (j!=0)*sin(omega2*j+delta2) + (j==0);
        AQ[j] = std::abs(Q[j]);
  	  }
    }

    auto L = P/Q;
    auto AL = AP/AQ;
    auto P1 = L*Q;
    // auto AP1 = AL*AQ; // UNUSED?!

    for(int j=0; j < 1<<order; ++j){
      assert(std::abs(P[j] - P1[j])/(std::abs(AL[j])+std::abs(AP[j])+std::abs(AQ[j])+std::numeric_limits<long double>::min()) < factor*std::numeric_limits<long double>::epsilon() );
    }
  }

  void test_division_as_multiplication(int order, long double factor,
				       long double omega, long double delta){

    NilpotentPolynomial<long double> P(order), Q(order), AP(order), AQ(order);
    for(int j=0; j < 1<<order; ++j){
      P[j] = sin(omega*j+delta);
      AP[j] = std::abs(P[j]);
    }

    //Q = P / [(1-(x_1+...+x_n))^{-1} ] = P * (1-(x_1+...+x_n))
    //(1-(x_1+...+x_n))^{-1} = 1 + (x_1+...+x_n) + ... + k! x_1 ... x_k +...+ n! x_1 ... x_n
    for(int j=0; j < 1<<order; ++j){
      Q[j] = ffd::core_math::Factorial(ffd::set_theory::CardinalitySet(j));
    }

    auto L1 = P/Q;

    NilpotentPolynomial<long double> Qm1(order);
    for(int j=0; j<order; ++j){

      Qm1[1<<j] = -1;
    }
    Qm1[0] = 1;

    auto L2 = P*Qm1;

    for(int j=0; j < 1<<order; ++j){
      assert(std::abs(L1[j] - L2[j])/std::abs(Q[j]) < factor*std::numeric_limits<long double>::epsilon() );
    }
  }

  void test_even_divide(int order, int factor, int omega1, int omega2, int delta1, int delta2){
	  NilpotentPolynomial<long double, Full> P(order), Q(order), AP(order), AQ(order);
		NilpotentPolynomial<long double, Even> P_(order), Q_(order), AP_(order), AQ_(order);
  	for(int j=0; j < 1<<order; ++j){
  	  if (__builtin_parity(j)==0){
  		  P[j] = sin(omega1*j+delta1);
  		  Q[j] = (j!=0)*sin(omega2*j+delta2) + (j==0);
        P_[j] = P[j];
        Q_[j] = Q[j];
  	  }
  	}

  	auto L = P/Q;
    auto L_ = P_/Q_;

  	for(int j=0; j < 1<<order; ++j){
  	  assert(std::abs(L[j] - L_[j]) < factor*std::numeric_limits<long double>::epsilon());
  	}
  }

  void operator_divide(){
    for(int k=0; k<12; ++k){

      long double factor = 10<<k;

      test_divide<>(k, factor, 3, 5, -12, -123);
      test_divide<>(k, factor, 13, 11, -532, -3283);
      test_divide<>(k, factor, 17, 19, -216, -60192);
      test_divide<>(k, factor, 23, 7, -816, -192);

  	  test_divide<Even>(k, factor, 3, 5, -12, -123);
  	  test_divide<Even>(k, factor, 13, 11, -532, -3283);
  	  test_divide<Even>(k, factor, 17, 19, -216, -60192);
  	  test_divide<Even>(k, factor, 23, 7, -816, -192);

      test_even_divide(k, factor, 3, 5, -12, -123);
      test_even_divide(k, factor, 13, 11, -532, -3283);
      test_even_divide(k, factor, 17, 19, -216, -60192);
      test_even_divide(k, factor, 23, 7, -816, -192);

      long double factor2 = 10;

      test_division_as_multiplication(k, factor2, 13, 2);
      test_division_as_multiplication(k, factor2, 17, -21);
      test_division_as_multiplication(k, factor2, 23, -121);
      test_division_as_multiplication(k, factor2, 17, -372);
      test_division_as_multiplication(k, factor2, 5, -823);
    }
  }


}//namespace
