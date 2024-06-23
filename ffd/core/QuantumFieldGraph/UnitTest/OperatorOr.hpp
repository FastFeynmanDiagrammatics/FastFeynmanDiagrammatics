

namespace ffd::qf_graph::unit_test{

  void OperatorOr(){
    int x = 3, y = 6;
    auto D1 = (Psi_(1)*Psi_(-1))(x);
    auto D2 = (Phi_(1)*Bar(Phi_(-1)))(y);
    auto G12 = D1|D2;

    static_assert(std::is_same<decltype(G12), QuantumFieldGraph>::value);
    assert(std::size(G12) == 2);

    G12 |= D1;
    QuantumFieldVertex V1(D1);
    auto GVD = V1|D2;
    auto GDV = D2|V1;

    static_assert(std::is_same<decltype(GVD), QuantumFieldGraph>::value);
    static_assert(std::is_same<decltype(GDV), QuantumFieldGraph>::value);
    assert(std::size(GVD) == 2);
    assert(std::size(GDV) == 2);

    GVD |= GDV;

    static_assert(std::is_same<decltype(GVD), QuantumFieldGraph>::value);
    assert(std::size(GVD) == 4);

    assert(GVD[1] == QuantumFieldVertex(D2));
    assert(GVD[2] == QuantumFieldVertex(D2));

    assert(GVD[3] == V1);
    assert(GVD[1][0] == D2);
    assert(GVD[3][0] == V1[0]);

    assert(std::any_cast<int>(GVD[3][0].Position) == x);
    assert(std::any_cast<int>(GVD[2][0].Position) == y);

    assert(GVD[3][0][0] == Psi_(1)[0]);
    assert(GVD[3][0][1] == Psi_(-1)[0]);
    assert(GVD[2][0][1] == Phi_(-1)[0]);
    
  }

}//namespace
