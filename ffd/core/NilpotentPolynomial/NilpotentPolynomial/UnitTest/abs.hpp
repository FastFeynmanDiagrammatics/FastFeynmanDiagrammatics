

namespace ffd::nilpotent_polynomial::unit_test{

  void absCheck(){
    NilpotentPolynomial<Real> P(5);
    for(BinaryInt S = 0; S < P.size(); ++S){
      P[S] = pow(-1.,S) * (S + 1.);
    }
    P = abs(P);
    for(BinaryInt S = 0; S < P.size(); ++S){
        assert( P[S] > 0.);
    }
  }
  
}//namespace
