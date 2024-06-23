

namespace ffd::qf_product::unit_test{

  void ProductCreation(){
    QuantumFieldProduct P_;
    P_.resize(10);

    assert(std::size(P_) == 10);

    QuantumField Q_;
    Q_.Dagger() = 1; Q_.Component() = 0; Q_.IsFermion() = true;
    QuantumFieldProduct P2_(Q_);

    assert(P2_[0] == Q_);
    
  }

}
