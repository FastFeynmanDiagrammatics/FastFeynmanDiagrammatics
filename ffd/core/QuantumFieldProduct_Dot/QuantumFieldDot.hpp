namespace ffd::qf_dot{

  using ffd::quantum_field::QuantumField;

  
  class QuantumFieldDot: public std::vector<QuantumField>{
  public:
    std::any Position;

    QuantumFieldDot(){}
    
  };

}//namespace ffd::qf_dot
