namespace ffd::user_space{
  
  ffd::quantum_field::QuantumField
  Nambu(ffd::quantum_field::QuantumField const& Q_){
    auto ret = Q_;
    if( ret.Component() < 0 ){
       ret.Dagger() *= -1;
    }
    ret.NotNambu() = !ret.NotNambu();
    return ret;
  }

}//namespace ffd
