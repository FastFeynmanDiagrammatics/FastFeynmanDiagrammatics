

namespace ffd::quantum_field::unit_test{

  void QuantumFieldCreation(){
    QuantumField G;
    G.Dagger() = 1;
    assert(G.Dagger() == 1);
  }

}//namespace
