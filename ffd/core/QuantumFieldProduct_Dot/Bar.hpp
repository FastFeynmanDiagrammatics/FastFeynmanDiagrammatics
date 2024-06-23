namespace ffd::user_space{

  // QuantumFieldProduct Bar(const QuantumFieldProduct& P_){
  //   QuantumFieldProduct ret = P_;
  //   for(int j=0; j < std::size(P_); ++j){
  //     ret[j] = Bar(P_[size(P_)-1-j]);
  //   }
  //   return ret;
  // }

  template<typename T>
  T Bar(T const& P_){
    T ret = P_;
    for(std::size_t j=0; j < size(P_); ++j){
      ret[j] = Bar(P_[size(P_)-1-j]);
    }
    return ret;
  }
  
}
