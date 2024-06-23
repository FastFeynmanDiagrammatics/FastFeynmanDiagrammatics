

namespace ffd::quantum_field::unit_test{

  void bar(){
    QuantumField Q1, Q2;
    Q1.Dagger() = 1;
    Q1.Component() = 1;
    Q1.IsFermion() = true;
    Q2 = Q1;
    Q1 = ffd::user_space::Bar(Q1);
    Q2.Dagger() *= -1;
    assert(Q1 == Q2);
  }

}
