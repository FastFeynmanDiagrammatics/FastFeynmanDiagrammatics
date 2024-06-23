namespace ffd::feynman_edge::unit_test{

  void quantum_field_position_function(){
    Real tau = 1.;
    auto dot = Psi_(1)(tau);
    auto qf_pos = quantum_field_position(dot);

    static_assert(std::is_same_v<decltype(qf_pos),
		  ffd::feynman_edge::QuantumFieldPosition>);



    Real tau2 = 0.;
    auto vertex = Psi_(1)(tau)*Bar(Psi_(-1))(tau2);
    auto qf_pos2 = quantum_field_positions(vertex);
    

    static_assert(std::is_same_v<decltype(qf_pos2),
		  ffd::feynman_edge::QuantumFieldPositions>);
    


  }

}//namespace
