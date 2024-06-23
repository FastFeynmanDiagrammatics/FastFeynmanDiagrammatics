

namespace ffd::user_space{

  using ffd::qf_vertex::QuantumFieldVertex;

  using ffd::qf_dot::QuantumFieldDot;
    
  QuantumFieldVertex operator*(const QuantumFieldVertex& V1_, const QuantumFieldVertex& V2_){
    QuantumFieldVertex ret;
    ret.resize( size(V1_) + size(V2_) );
    for(std::size_t j=0; j < size(V1_); ++j){
      ret[j] = V1_[j];
    }
    for(std::size_t j=0; j < size(V2_); ++j){
      ret[j + size(V1_)] = V2_[j];
    }
    return ret;
  }

}//namespace ffd::user_space
