

namespace ffd::user_space{

  using ffd::qf_graph::QuantumFieldGraph;
  
  QuantumFieldGraph operator|(const QuantumFieldGraph& G1_, const QuantumFieldGraph& G2_){
    QuantumFieldGraph ret;
    ret.resize( size(G1_) + size(G2_) );
    for(std::size_t j=0; j < size(G1_); ++j){
      ret[j] = G1_[j];
    }
    for(std::size_t j=0; j < size(G2_); ++j){
      ret[j+size(G1_)] = G2_[j];
    }
    return ret;
  }
  
}//namespace ffd::user_space


namespace ffd::qf_graph{

  QuantumFieldGraph& QuantumFieldGraph::operator|=(const QuantumFieldGraph& G_){
    using ffd::user_space::operator|;
    return *this = *this | G_;
  }
  
}//namespace ffd::qf_graph
