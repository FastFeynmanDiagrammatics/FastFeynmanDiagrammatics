namespace ffd::user_space{

  using ffd::quantum_field::QuantumField;

  using ffd::qf_product::QuantumFieldProduct;
  
  QuantumFieldProduct Psi_(int Component_){
    QuantumField ret;
    ret.IsFermion() = true;
    ret.Dagger() = 1;
    ret.Component() = Component_;
    return ret;
  }


  
  QuantumFieldProduct NambuPsi_(int Component_){
    auto ret = Psi_(Component_);
    ret[0].NotNambu() = false;
    if( Component_ < 0 ){
      ret[0].Dagger() = -1;
    }
    return ret;
  }

  
  QuantumFieldProduct Psi_(const char* Spin_){
    if(std::string(Spin_) == std::string("up")){
      return Psi_(1);
    }else if(
	     std::string(Spin_) == std::string("down") ||
	     std::string(Spin_) == std::string("do") ||
	     std::string(Spin_) == std::string("dn")
	     ){
      return Psi_(-1);
    }else{
      assert(false);
      return Psi_(0);
    }
  }

  
  QuantumFieldProduct Rho_(int Component_){
    QuantumField ret;
    ret.IsFermion() = true;
    ret.Dagger() = 0;
    ret.Component() = Component_;
    return ret;
  }


  QuantumFieldProduct Phi_(int Component_){
    QuantumField ret;
    ret.IsFermion() = false;
    ret.Dagger() = 0;
    ret.Component() = Component_;
    return ret;
  }


  QuantumFieldProduct Eta_(int Component_){
    QuantumField ret;
    ret.IsFermion() = false;
    ret.Dagger() = 1;
    ret.Component() = Component_;
    return ret;
  }


}//namespace ffd::quantum_field
