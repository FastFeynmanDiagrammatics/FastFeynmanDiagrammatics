

namespace ffd::qf_graph::unit_test{

  void IteratorCheck(){
    int x = 3;
    auto D1 = (Psi_(1)*Psi_(-1))(x);
    auto D2 = (Bar(Psi_(1)*Psi_(-1)))(x)*D1;
    auto D3 = (Bar(Psi_(1)*Psi_(-1))*FlipSpin(Phi_(0)))(x);
    auto D4 = (Rho_(3)*Phi_(-1)*Psi_(3)*Psi_(1)*Bar(Phi_(2)))(x);
    auto G = D1|D2;
    G |= D3|D4;
    QFGraphIterator it(G);
    auto res = it();

    //checking operator()
    auto& [QF, pos] = res;
    assert(QF == Psi_(1)[0]);
    auto& [vertex, dot] = pos;
    assert(vertex == 0);
    assert(dot == 0);

    
    //checking operator++
    res = (++it)();
    assert(QF == Psi_(-1)[0]);
    assert(vertex == 0);
    assert(dot == 0);

    res = (++it)();
    assert(QF == Bar(Psi_(-1))[0]);
    assert(vertex == 1);
    assert(dot == 0);

    res = (++it)();
    assert(QF == Bar(Psi_(1))[0]);
    assert(vertex == 1);
    assert(dot == 0);

    res = (++it)();
    assert(QF == Psi_(1)[0]);
    assert(vertex == 1);
    assert(dot == 1);

    res = (++it)();
    assert(QF == Psi_(-1)[0]);
    assert(vertex == 1);
    assert(dot == 1);

    res = (++it)();
    assert(QF == Bar(Psi_(-1))[0]);
    assert(vertex == 2);
    assert(dot == 0);

    res = (++it)();
    assert(QF == Bar(Psi_(1))[0]);
    assert(vertex == 2);
    assert(dot == 0);

    QFGraphIterator it2(G);
    it = it2 + 1;
    res = it();
    assert(QF == Psi_(-1)[0]);
    assert(vertex == 0);
    assert(dot == 0);

    //checking NotEnd()
    it2 = it2 + 13;
    assert(it2.NotEnd());
    ++it2;
    assert(!it2.NotEnd());
    
  }

}//namespace
