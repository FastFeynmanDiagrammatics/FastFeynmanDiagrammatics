

namespace ffd::qf_product::unit_test{

  void ProductOperatorPar(){
    QuantumField Q_; Q_.Dagger()=1; Q_.Component()=1; Q_.IsFermion()=false;
    QuantumFieldProduct P_(Q_);
    P_.push_back(ffd::user_space::Bar(Q_));
    P_.push_back(ffd::user_space::FlipSpin(Q_));

    auto D_ = P_(static_cast<int>(111));

    assert(std::any_cast<int>(D_.Position) == 111);

    for(std::size_t j=0; j< std::size(P_); ++j){
      assert(P_[j] == D_[j]);
    }
    
  }

}
