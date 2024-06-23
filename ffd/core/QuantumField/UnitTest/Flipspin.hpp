

namespace ffd::quantum_field::unit_test{

  void Flipspin(){
    QuantumField Q1, Q2;
    Q1.Dagger() = 1;
    Q1.Component() = 1;
    Q1.IsFermion() = true;
    Q2 = Q1;
    Q1 = ffd::user_space::FlipSpin(Q1);
    Q2.Component() *= -1;
    assert(Q1 == Q2);
  }

}
