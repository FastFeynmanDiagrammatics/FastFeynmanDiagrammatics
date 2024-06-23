

namespace ffd::user_space{

  // QuantumFieldProduct FlipSpin(const QuantumFieldProduct& Q_){
  //   QuantumFieldProduct ret = Q_;
  //   for(int j=0; j < std::size(ret); ++j){
  //     ret[j].Component() *= -1;
  //   }
  //   return ret;
  // }

  template<typename T>
  T FlipSpin(T const& X_){
    T ret = X_;
    for(std::size_t j=0; j < std::size(ret); ++j){
      ret[j] = FlipSpin(ret[j]);
    }
    return ret;
  }
  
}//namespace ffd::user_space
