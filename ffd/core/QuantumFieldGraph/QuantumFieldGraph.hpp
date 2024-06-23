namespace ffd::qf_graph{

  using ffd::quantum_field::QuantumField;

  using ffd::qf_dot::QuantumFieldDot;
  
  using ffd::qf_vertex::QuantumFieldVertex;

  
  class QuantumFieldGraph: public std::vector<QuantumFieldVertex>{
  public:
    
    QuantumFieldGraph(){}

    QuantumFieldGraph(QuantumFieldVertex const& V_){
      this->push_back(V_);
    }
    
    QuantumFieldGraph(const QuantumFieldDot& D_): QuantumFieldGraph(QuantumFieldVertex(D_)) {}
    
    QuantumFieldGraph& operator|=(const QuantumFieldGraph&);

    //used to feed WickFunction
    // std::pair<QuantumFieldVertex const&, char> operator()(char j_){
    //   return {std::ref((*this)[j_]), j_};
    // }
    std::pair<QuantumFieldVertex, char> operator()(char j_){
      return {(*this)[j_], j_};
    }

  };

}//namespace ffd::qf_graph
