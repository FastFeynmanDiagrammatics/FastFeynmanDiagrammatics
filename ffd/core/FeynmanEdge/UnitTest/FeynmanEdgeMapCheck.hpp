

namespace ffd::feynman_edge::unit_test{

  using namespace ffd::user_space;
  
  void FeynmanEdgeMapCheck(){
    int x = 1, y = 0;
    auto D1 = (Psi_(1)*Psi_(-1)*Phi_(0))(x);
    auto D2 = (Bar(Psi_(1)*Psi_(-1)*Phi_(0)))(y);
    auto G = D1|D2;

    auto drawer = [](QuantumFieldPositions Q)->int{
		    auto [x_inc, x_out] = Q;
		    return -std::any_cast<int>(x_inc.second)+std::any_cast<int>(x_out.second);};

    auto map = CreateFeynmanEdgeMap<int>(G, drawer);

    using std::size;
    assert(size(map) == 3);

    ffd::qf_graph::QFGraphIterator it(G);
    
    assert(map.count({(it+5)(), it()}) == 1);
    assert(map.count({(it)(), (it+5)()}) == 0);
    assert(map.count({(it+4)(), (it+1)()}) == 1);
    assert(map.count({(it+2)(), (it+3)()}) == 1);

    assert(map.at({(it+5)(), it()}) == 1);
    assert(map.at({(it+4)(), (it+1)()}) == 1);
    assert(map.at({(it+2)(), (it+3)()}) == -1);

  }

}//namespace 
