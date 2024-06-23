namespace ffd::qf_dot{

  auto
  quantum_field_position(ffd::qf_dot::QuantumFieldDot const& dot_){
    assert( size(dot_) == 1 );
    std::any X = dot_.Position;

    
    return std::make_pair(dot_[0], X);
  }
  
}//namespace ffd::qf_dot



namespace ffd::qf_vertex{
  
  auto
  quantum_field_positions(ffd::qf_vertex::QuantumFieldVertex const& vertex_){
    assert( size(vertex_) == 2 );
    auto ar_0 = quantum_field_position(vertex_[0]);
    auto ar_1 = quantum_field_position(vertex_[1]);
    std::array<decltype(ar_0), 2> ret{ar_0, ar_1};

    
    return ret;
  }

}//namespace ffd::qf_vertex{
