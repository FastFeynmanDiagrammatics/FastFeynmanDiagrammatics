

namespace ffd::qf_graph::unit_test{

  void OperatorPar(){
    float x = 1.;
    auto D = (Psi_(1)*Psi_(-1))(x);
    QuantumFieldVertex V(D);
    QuantumFieldGraph G(D);
    G.push_back(D);

    assert(G(0).first == V);
    assert(G(0).second == 0);
    assert((G(1) == std::pair{V, static_cast<decltype(G(1).second)>(1)}));
    
  }

}
