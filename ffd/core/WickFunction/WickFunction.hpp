namespace ffd::wick_function{

  using ffd::quantum_field::QuantumField;
  using ffd::qf_dot::QuantumFieldDot;
  using ffd::feynman_edge::QuantumFieldVertexDot;
  using ffd::qf_vertex::QuantumFieldVertex;
  using ffd::qf_graph::QuantumFieldGraph;

  
  
  class WickFunction: public
  std::vector<std::array<std::vector<ffd::feynman_edge::QuantumFieldVertexDot>, 2>>{
  public:

    WickFunction(std::function<char(QuantumField)> WhichBlock_ =
		 ffd::feynman_edge::WhichBlockDefault):
      WhichBlock(WhichBlock_) {}
    
    std::vector<char> BlockNumbers;

    std::vector<bool> IsFermionicBlock;

    std::function<char(QuantumField)> WhichBlock;
    
    char Sign = 1;

    void reset(){
      this->resize(0);
      BlockNumbers.resize(0);
      IsFermionicBlock.resize(0);
      Sign = 1;
    }

    WickFunction& operator*=(std::pair<QuantumFieldVertex, char>);
    
  };


  
}//namespace ffd::wick_function
