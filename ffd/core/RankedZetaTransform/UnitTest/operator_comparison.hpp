namespace ffd::ranked_zeta::unit_test{

	void comparison(){
		using NilpotentPolynomial = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Full>;
		using namespace user_space;

      const int order = 6;
      NilpotentPolynomial P(order);
      for(int j=0; j < 1<<order; ++j){
        P[j] = j;
      }
	  RankedZeta Z = P;
      auto R = Z;

	  double x = 0.0;

	  assert(R==Z);
	  assert(R>=Z);
	  assert(R<=Z);
	  assert(x>=Z);
	  assert(x<=Z);
	  assert(Z>=x);
	  assert(Z<=x);

      x = 1.23456789238123876213l;
      R += x;

	  assert(R>Z);
	  assert(Z<R);
	  assert(x>Z);
	  assert(Z<x);

      x = -1.23456789238123876213l;

	  assert(x<Z);
	  assert(Z>x);

    }

}//namespace
