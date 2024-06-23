namespace ffd::nilpotent_polynomial::unit_test{


  void test_multiply_order4(long double omega1, long double omega2,
			    long double delta1, long double delta2, long double factor = 2){
    const int order = 4;
    NilpotentPolynomial<long double> P(order), Q(order), AP(order), AQ(order);
    for(int j=0; j < 1<<order; ++j){
      P[j] = std::sin(omega1*j+delta1);
      Q[j] = std::sin(omega2*j+delta2);
      AP[j] = std::abs(P[j]);
      AQ[j] = std::abs(Q[j]);
    }

    NilpotentPolynomial<long double> L;
    L *= P;

    assert(!L.NotZero);
    using std::size;
    assert(size(L) == 1);

    const long double x = 1.23456789123123131231231242l;
    L = x;
    assert(L.NotZero);
    assert(size(L) == 1);

    L *= P;

    assert(L.NotZero);
    assert(size(L) == 1<<order);
    for(int j=0; j < 1<<order; ++j){
      assert(std::abs(L[j] - P[j]*x) < std::abs(L[j])*std::numeric_limits<long double>::epsilon() );
    }

    L = P*Q;
    auto AL = AP*AQ;

    assert(std::abs(L[0] - P[0]*Q[0]) < std::abs(L[0])*std::numeric_limits<long double>::epsilon() );
    assert(std::abs(L[1] - P[1]*Q[0] - P[0]*Q[1]) < std::abs(L[1])*std::numeric_limits<long double>::epsilon() );
    assert(std::abs(L[0b10] - P[0b10]*Q[0] - P[0]*Q[0b10])
	   < std::abs(AL[0b10])*std::numeric_limits<long double>::epsilon() );
    assert(std::abs(L[0b100] - P[0b100]*Q[0] - P[0]*Q[0b100])
	   < std::abs(AL[0b100])*std::numeric_limits<long double>::epsilon() );
    auto L0011 = P[0b11]*Q[0] + P[0b10]*Q[0b01] + P[0b01]*Q[0b10] + P[0]*Q[0b11];
    assert(std::abs(L[0b11] - L0011) < std::abs(AL[0b0011])*std::numeric_limits<long double>::epsilon() );
    auto L0110 = P[0b110]*Q[0] + P[0b100]*Q[0b010] + P[0b010]*Q[0b100] + P[0]*Q[0b110];
    assert(std::abs(L[0b110] - L0110) < std::abs(AL[0b0110])*std::numeric_limits<long double>::epsilon() );
    auto L0101 = P[0b101]*Q[0] + P[0b100]*Q[0b001] + P[0b001]*Q[0b100] + P[0]*Q[0b101];
    assert(std::abs(L[0b101] - L0101) < std::abs(AL[0b0101])*std::numeric_limits<long double>::epsilon() );
    auto L0111 = P[0b111]*Q[0] + P[0b110]*Q[0b001] + P[0b101]*Q[0b010] + P[0b011]*Q[0b100]+
      P[0b100]*Q[0b011] + P[0b010]*Q[0b101] + P[0b001]*Q[0b110] + P[0]*Q[0b111];
    assert(std::abs(L[0b111] - L0111) < std::abs(AL[0b0111])*std::numeric_limits<long double>::epsilon() );
    auto L1011 = P[0b1011]*Q[0] + P[0b1010]*Q[0b0001] + P[0b1001]*Q[0b0010] + P[0b0011]*Q[0b1000]+
      P[0b1000]*Q[0b0011] + P[0b0010]*Q[0b1001] + P[0b0001]*Q[0b1010] + P[0]*Q[0b1011];
    assert(std::abs(L[0b1011] - L1011) < std::abs(AL[0b1011])*std::numeric_limits<long double>::epsilon() );
    auto L1101 = P[0b1101]*Q[0] + P[0b1100]*Q[0b001] + P[0b1001]*Q[0b0100] + P[0b0101]*Q[0b1000]+
      P[0b1000]*Q[0b0101] + P[0b0100]*Q[0b1001] + P[0b0001]*Q[0b1100] + P[0]*Q[0b1101];
    assert(std::abs(L[0b1101] - L1101) < std::abs(AL[0b1101])*std::numeric_limits<long double>::epsilon() );
    auto L1110 = P[0b1110]*Q[0] + P[0b1100]*Q[0b0010] + P[0b1010]*Q[0b0100] + P[0b0110]*Q[0b1000]+
      P[0b1000]*Q[0b0110] + P[0b0100]*Q[0b1010] + P[0b0010]*Q[0b1100] + P[0]*Q[0b1110];
    assert(std::abs(L[0b1110] - L1110) < std::abs(AL[0b1110])*std::numeric_limits<long double>::epsilon() );
    auto L1111 = P[0b1111]*Q[0] + P[0b1110]*Q[0b0001] + P[0b1101]*Q[0b0010] + P[0b1011]*Q[0b0100]+
      P[0b0111]*Q[0b1000] + P[0b1100]*Q[0b0011] + P[0b1010]*Q[0b0101] + P[0b1001]*Q[0b0110] +
      P[0b0110]*Q[0b1001] + P[0b0101]*Q[0b1010] + P[0b0011]*Q[0b1100] + P[0b1000]*Q[0b0111] +
      P[0b0100]*Q[0b1011] + P[0b0010]*Q[0b1101] + P[0b0001]*Q[0b1110] + P[0]*Q[0b1111];
    assert(std::abs(L[0b1111] - L1111) < factor*std::abs(AL[0b1111])*std::numeric_limits<long double>::epsilon() );

  }

  void test_even_multiply(int order, int factor, int omega1, int omega2, int delta1, int delta2){
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

  	auto L = P*Q;
    auto L_ = P_*Q_;

  	for(int j=0; j < 1<<order; ++j){
  	  assert(std::abs(L[j] - L_[j]) < factor*std::numeric_limits<long double>::epsilon());
  	}
  }


  void operator_multiply(){
    test_multiply_order4(13, 5, 0.2, 0.4);
    test_multiply_order4(13, 7, 0.33, 0.25);
    test_multiply_order4(17, 9, 0.82123, 0.123);
    test_multiply_order4(11, 3, 0.2123, 0.7777);
    test_multiply_order4(23, 19, 2.123, -4.212);
    test_multiply_order4(29, 37, 8.921, 5.1231);
    test_multiply_order4(2, 43, -321., -1237);
    test_multiply_order4(3, 31, -1234, -98732);
    test_multiply_order4(3, 3, -1234, -1234, 2);  //does not pass with factor == 1
    test_multiply_order4(2, 2, 1, 1);
    test_multiply_order4(0, 0, 1, 1);
    test_multiply_order4(0, 1, 1, 1);

	for(int k=0; k<12; ++k){

	  long double factor = 10<<k;
		test_even_multiply(k, factor, 3, 5, -12, -123);
		test_even_multiply(k, factor, 13, 11, -532, -3283);
		test_even_multiply(k, factor, 17, 19, -216, -60192);
		test_even_multiply(k, factor, 23, 7, -816, -192);
	}
  }

}//namespace
