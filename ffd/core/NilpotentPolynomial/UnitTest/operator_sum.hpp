

namespace ffd::nilpotent_polynomial::unit_test{

  void operator_sum(){
    const int order = 10;
    NilpotentPolynomial P(order);
    for(int j=0; j < 1<<order; ++j){
      P[j] = 1 + std::log(1+j);
    }
    static_assert(std::is_same<decltype(P), NilpotentPolynomial<Real>>::value);

    NilpotentPolynomial<Real> Q;
    Q = P + Q;
    
    for(int j=0; j < 1<<order; ++j){
      assert(std::abs(P[j]-Q[j]) < 10*std::numeric_limits<Real>::epsilon());
    }

    auto L = P - Q;

    for(int j=0; j < 1<<order; ++j){
      assert(std::abs(L[j]) < 10*std::numeric_limits<Real>::epsilon());
    }
    

  }

}//namespace
