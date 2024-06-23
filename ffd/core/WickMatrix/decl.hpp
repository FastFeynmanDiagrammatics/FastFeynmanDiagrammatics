namespace ffd::wick_matrix{

  using ffd::quantum_field::QuantumField;

  using ffd::qf_product::QuantumFieldProduct;

  using intTypeDefault = int;

  
  template<typename Field, typename intType = intTypeDefault>
  class WickMatrix{
  public:
    std::vector<Field> Components;

    intType ConvertToInternalRepresentation(intType, intType) const;

    intType NumberDestructionOperators = 0;

    bool IsFermion = true;
    bool NotHermitian = true;

    WickMatrix(){};

    WickMatrix(std::pair<bool, bool>, intType);

    WickMatrix(QuantumField const&, intType);

    WickMatrix(QuantumFieldProduct const&, intType);

    Field operator()(intType, intType) const;
    Field& operator()(intType, intType, const char);

    template<typename T1, typename T2>
    friend std::size_t size(WickMatrix<T1, T2> const&);
    
  };

  
  
  template<typename Field, typename intType>
  std::size_t size(WickMatrix<Field, intType> const& M_){
    return M_.NumberDestructionOperators;
  }

  
}//namespace ffd::wick_matrix
