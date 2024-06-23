namespace ffd::nilpotent_polynomial::unit_test{

	void comparison(){

      const int order = 6;
      NilpotentPolynomial<Real> P(order);
      for(int j=0; j < 1<<order; ++j){
        P[j] = j;
      }
      auto Q = P;

	  double x = 0.0;

	  assert(Q==P);
	  assert(Q>=P);
	  assert(Q<=P);
	  assert(x>=P);
	  assert(x<=P);
	  assert(P>=x);
	  assert(P<=x);

      x = 1.23456789238123876213l;
      Q += x;

	  assert(Q>P);
	  assert(P<Q);
	  assert(x>P);
	  assert(P<x);

      x = -1.23456789238123876213l;

	  assert(x<P);
	  assert(P>x);

    }

}//namespace
