namespace ffd::user_space{

  std::ostream& operator<<(std::ostream& out, const ffd::quantum_field::QuantumField& Q_){
    if(Q_.Dagger() == -1){
      out<<"Bar(";
    }
    if(Q_.NotNambu() == false){
      out<<"Nambu";
    }
    if(Q_.IsFermion()){
      if(Q_.Dagger() != 0){
	out<<"Psi_";
      }else{
	out<<"Rho_";
      }
    }else{
      if(Q_.Dagger() != 0){
	out<<"Eta_";
      }else{
	out<<"Phi_";
      }
    }
    out<<"{"<<int(Q_.Component())<<"}";
    if(Q_.Dagger() == -1){
      out<<")";
    }
    return out;
  }
  

}//namespace ffd::quantum_field
