

namespace ffd::user_space{

  using ffd::qf_product::QuantumFieldProduct;
  
  QuantumFieldProduct operator*(const QuantumFieldProduct& P1_, const QuantumFieldProduct& P2_){
    QuantumFieldProduct ret;
    ret.resize(size(P1_) + size(P2_));
    for(std::size_t j=0; j < size(P1_); ++j){
      ret[j] = P1_[j];
    }
    for(std::size_t j=0; j < size(P2_); ++j){
      ret[j + size(P1_)] = P2_[j];
    }
    return ret;
  }
  

  
  QuantumFieldProduct&
  operator*=(QuantumFieldProduct& P1_, const QuantumFieldProduct& P2_){
    P1_ = P1_*P2_;
    return P1_;
  }
  

  
}
