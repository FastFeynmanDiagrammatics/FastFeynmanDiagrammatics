namespace ffd::user_space{

  template<typename Field>
  class QuadraticActionTerm: public std::array<ffd::qf_dot::QuantumFieldDot, 2>{
  public:
    Field Value;
    
    QuadraticActionTerm(){}
    
  };

}//namespace ffd::user_space
