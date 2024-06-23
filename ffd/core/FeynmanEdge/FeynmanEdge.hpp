namespace ffd::feynman_edge{

  using ffd::quantum_field::QuantumField;

  using QuantumFieldVertexDot = std::pair<QuantumField, std::array<char, 2>>;

  using QuantumFieldPosition  = std::pair<QuantumField, std::any>;
  
  using FeynmanEdge           = std::array<QuantumFieldVertexDot, 2>;

  using QuantumFieldPositions = std::array<QuantumFieldPosition, 2>;
  
  template<typename Field>
   // using FeynmanEdgeMap = std::map<FeynmanEdge, Field>;
  using FeynmanEdgeMap = ffd::flat_map::FlatMap<FeynmanEdge, Field>;
  
}//namespace ffd::feynman_edge

