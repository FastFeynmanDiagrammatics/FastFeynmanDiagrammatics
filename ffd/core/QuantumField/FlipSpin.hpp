

namespace ffd::user_space{

  auto FlipSpin(ffd::quantum_field::QuantumField const& Q_){
    auto ret = Q_;
    ret.Component() *= -1;
    return ret;
  }

}//namespace ffd::user_space
