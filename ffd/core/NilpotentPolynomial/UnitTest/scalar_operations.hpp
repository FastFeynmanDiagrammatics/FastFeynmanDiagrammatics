

namespace ffd::nilpotent_polynomial::unit_test{

  void scalar_operations(){
    
    const int order = 10;
    NilpotentPolynomial<Real> P(order);
    for(int j=0; j < 1<<order; ++j){
      P[j] = j;
    }

    auto Q = P;
    long double x = 1.23456789238123876213l;
    Q += x;

    for(int j=0; j < 1<<order; ++j){
      assert(std::abs(Q[j] - (j + x*(j==0))) < 10*std::numeric_limits<Real>::epsilon());
    }

    Q = P - x;

    for(int j=0; j < 1<<order; ++j){
      assert(std::abs(Q[j] - (P[j] - x*(j==0))) < 10*std::numeric_limits<Real>::epsilon());
    }

    Q = P;
    Q *= x;

    for(int j=0; j < 1<<order; ++j){
      assert(std::abs(Q[j] - x*P[j]) < 10*(std::numeric_limits<Real>::min()+std::abs(Q[j])*std::numeric_limits<Real>::epsilon()));
    }

    Q = P/x;

    for(int j=0; j < 1<<order; ++j){
      assert(std::abs(Q[j] - P[j]/x) < 10*(std::numeric_limits<Real>::min()+std::abs(Q[j])*std::numeric_limits<Real>::epsilon()));
    }
    
    Q = -P;
    
    for(int j=0; j < 1<<order; ++j){
      assert(std::abs(Q[j] + P[j]) < 10*std::numeric_limits<Real>::epsilon());
    }

  }

}//namespace 
