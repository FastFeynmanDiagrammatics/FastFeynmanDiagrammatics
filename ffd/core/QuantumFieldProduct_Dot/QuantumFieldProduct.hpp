namespace ffd::qf_product{

  using ffd::quantum_field::QuantumField;

  using ffd::qf_dot::QuantumFieldDot;

  
  class QuantumFieldProduct: public std::vector<QuantumField>{
  public:
    
    QuantumFieldProduct(){}
    
    QuantumFieldProduct(const QuantumField& Q_){
      this->resize(1);
      this->operator[](0) = Q_;
    }

    template<typename T>
    QuantumFieldDot operator()(T x){
      QuantumFieldDot ret;
      using std::size;
      ret.resize(size(*this));
      ret.Position = x;
      for(std::size_t j=0; j < size(*this); ++j){
	ret[j] = this->operator[](j);
      }
      return ret;
    }
    
  };

}
