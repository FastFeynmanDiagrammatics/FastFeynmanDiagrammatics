namespace ffd::nilpotent_polynomial::unit_test{

  void Creation(){
    NilpotentPolynomial<Real> P;

    assert(!P.NotZero);
    assert(P.size() == 1);
    assert(std::abs(P[0]) < 10*std::numeric_limits<Real>::min() );
    //    std::cerr << "here";
    const int order = 5;
    NilpotentPolynomial<float> P2(order);

    assert(!P2.NotZero);
    P2[0] = 0.;
    assert(P2.NotZero);

    assert(P2.size() == 1<<order);
    assert(P2.linear_size() == order);

    for(int j=0; j < 1<<order; ++j){
      assert(std::abs(P2[j]) < 10*std::numeric_limits<Real>::min() );
    }

    const long double z = 3.12312312312392438232323423748237l;
    NilpotentPolynomial<Real> P3 = z;
    assert(P3.NotZero);
    assert(P3.size() == 1);
    assert(std::abs(P3[0]-z) < 10*std::numeric_limits<Real>::epsilon() );


    std::vector<NilpotentPolynomial<Real>> P_vec;
  }

}
