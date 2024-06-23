

namespace ffd::quantum_field::unit_test{

  void OperatorLess(){
    QuantumField Q1, Q2;
    Q1.Dagger() = 1;
    Q2.Dagger() = 1;
    Q1.Component() = 1;
    Q2.Component() = -1;
    Q1.IsFermion() = false;
    Q2.IsFermion() = true;
    
    assert(Q2 < Q1);
  }

}
