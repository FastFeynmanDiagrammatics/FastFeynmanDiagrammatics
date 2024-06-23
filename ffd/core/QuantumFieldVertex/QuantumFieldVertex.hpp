namespace ffd::qf_vertex{

  
  using ffd::qf_dot::QuantumFieldDot;

  
  class QuantumFieldVertex: public std::vector<QuantumFieldDot>{
  public:

    std::optional<std::variant<char, std::string>> Type;
    
    QuantumFieldVertex() {}
    
    QuantumFieldVertex(QuantumFieldDot const& D_){
      (*this).push_back(D_);
    }
    
  };

  

}//namespace ffd::qf_vertex
