namespace ffd::user_space{

  using ffd::feynman_edge::QuantumFieldPositions;
  
  //we create the FeynmanEdgeMap starting from a QuantumFieldGraph,
  //a drawer, and a function that tells us to which block the QuantumField
  //belongs to.
  template<typename Field,
	   typename function_qfposition_field>
  ffd::feynman_edge::FeynmanEdgeMap<Field>
  CreateFeynmanEdgeMap(ffd::qf_graph::QuantumFieldGraph const& G_,
		       function_qfposition_field drawer_,
		       // std::function<Field(QuantumFieldPositions)> const& drawer_,
		       std::function<char(QuantumField)> WhichBlock_ =
		       ffd::feynman_edge::WhichBlockDefault,
		       bool NoSelfInteraction = true){
    using namespace ffd::feynman_edge;
    using namespace ffd::qf_graph;
    FeynmanEdgeMap<Field> ret;
    for(QFGraphIterator It1(G_); It1.NotEnd(); ++It1){
      for(QFGraphIterator It2 = It1+1; It2.NotEnd(); ++It2){
	auto const& QuantumFieldVertexDotMin  = It1() < It2() ? It1() : It2();
	auto const& [QuantumFieldMin, VertexDotMin] = QuantumFieldVertexDotMin;
	auto const& QuantumFieldVertexDotMax  = It1() < It2() ? It2() : It1();
	auto const& [QuantumFieldMax, VertexDotMax] = QuantumFieldVertexDotMax;
	if(WhichBlock_(QuantumFieldMin) == WhichBlock_(QuantumFieldMax)
	   && QuantumFieldMin.Dagger() == - QuantumFieldMax.Dagger()){
	  if(It1.Vertex != It2.Vertex || !NoSelfInteraction){
	    FeynmanEdge Edge = {QuantumFieldVertexDotMin, QuantumFieldVertexDotMax};
	    auto const& [VertexMin, DotMin] = VertexDotMin;
	    std::any PositionMin = G_[VertexMin][DotMin].Position;
	    auto const& [VertexMax, DotMax] = VertexDotMax;
	    std::any PositionMax = G_[VertexMax][DotMax].Position;
	    if(ret.count(Edge) == 0 ){
	      std::pair<QuantumField, std::any> const&
		QuantumFieldPositionMin = {QuantumFieldMin, PositionMin};
	      std::pair<QuantumField, std::any> const&
		QuantumFieldPositionMax = {QuantumFieldMax, PositionMax};
	      ret.insert_or_assign(Edge, drawer_({QuantumFieldPositionMin, QuantumFieldPositionMax}));
	    }
	  }else{
	    FeynmanEdge Edge = {QuantumFieldVertexDotMin, QuantumFieldVertexDotMax};
	    if(ret.count(Edge) == 0){
	      ret.insert_or_assign(Edge, 0.);
	    }
	  }
	}
      }
    }
    return ret;
  }

}//namespace ffd::qf_graph
