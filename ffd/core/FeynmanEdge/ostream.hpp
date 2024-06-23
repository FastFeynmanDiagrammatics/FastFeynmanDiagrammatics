namespace ffd::user_space{

  std::ostream& operator<<(std::ostream& out_, ffd::feynman_edge::QuantumFieldVertexDot QFD_){
    auto [QField, VD] = QFD_;
    auto [Vertex, Dot] = VD;
    return out_<<QField<<"("<<int(Vertex)<<", "<<int(Dot)<<")";
  }
  
  
  std::ostream& operator<<(std::ostream& out_, const ffd::feynman_edge::FeynmanEdge& E_){
    auto [x1, x2] = E_;
    return out_<<x1<<" -> "<<x2;
  }

  
  template<typename Field, typename std::enable_if_t<
			     std::is_same<ffd::feynman_edge::FeynmanEdgeMap<Field>,
					  std::map<ffd::feynman_edge::FeynmanEdge, Field>>::value, int > = 0>
  std::ostream& operator<<(std::ostream& out_, const ffd::feynman_edge::FeynmanEdgeMap<Field>& Map_){
    out_<<SetOstreamFloatFlags();
    for(auto [edge, value]: Map_){
      out_<<edge<<" = "<<value<<"\n";
    }
    return out_;
  }

}
