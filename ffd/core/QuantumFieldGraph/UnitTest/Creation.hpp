

namespace ffd::qf_graph::unit_test{

  using namespace ffd::user_space;
  
  void Creation(){
    QuantumFieldGraph G1;
    
    double x = 12, y= -12;
    auto V_ = (Psi_(1)*Psi_(-1)*Bar(Eta_(2)))(x)*Psi_(1)(y);
    QuantumFieldGraph G2(V_);
    
    auto D_ = (Psi_(1)*Psi_(-1)*Bar(Eta_(2)))(x);
    QuantumFieldGraph G3(D_);

    G1.push_back(V_);
    G1.push_back(V_);

    assert(G2[0] == V_);
    assert(G1[1] == V_);
    assert(std::size(G1) == 2);
    assert(std::size(G3) == 1);

    QuantumFieldVertex VD_(D_);

    assert(G3[0] == VD_);
    
  }

}
