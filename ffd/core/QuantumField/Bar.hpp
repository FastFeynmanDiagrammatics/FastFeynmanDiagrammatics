namespace ffd::user_space{
  
  ffd::quantum_field::QuantumField
  Bar(ffd::quantum_field::QuantumField Q_){
    auto ret = Q_;
    ret.Dagger() *= -1;
    return ret;
  }

}//namespace ffd


