

namespace ffd::qf_vertex::unit_test{

  using namespace ffd::user_space;
  
  void Multiply(){
    double x = 11;
    int j = 12;
    auto Q_ = (Psi_(1)*FlipSpin(Psi_(-1)))(x)*(Eta_(0)*Phi_(-2)*Bar(Psi_(2)))(j);

    static_assert(std::is_same<decltype(Q_), QuantumFieldVertex>::value);
    
    static_assert(std::is_same<decltype(std::any_cast<double>(Q_[0].Position)), double>::value);
    static_assert(std::is_same<decltype(std::any_cast<int>(Q_[1].Position)), int>::value);

    assert(std::size(Q_) == 2);
    assert(std::size(Q_[0]) == 2);
    assert(std::size(Q_[1]) == 3);
    assert(Q_[0][1].Component() == 1);
    assert(Q_[1][2].Dagger() == -1);

    auto Q2_ = Q_*Q_*Q_;

    assert(std::size(Q2_) == 6);
    assert(Q2_[1] == Q2_[3]);
    assert(Q2_[1] == Q2_[5]);

  }

  
}//namespace
