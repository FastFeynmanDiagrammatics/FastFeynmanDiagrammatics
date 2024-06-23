namespace ffd::core_integration_test{

  using ffd::feynman_edge::QuantumFieldPositions;

  template<typename Field>
  class BarePropagator{
  private:
    Field beta, mu;
  public:
    BarePropagator(Field InverseTemperature, Field ChemicalPotential):
      beta(InverseTemperature), mu(ChemicalPotential){}

    Field operator()(QuantumFieldPositions X) const{
      auto const& [XMin, XMax] = X;
      auto const& [QuantumFieldMin, PositionMin] = XMin;
      auto const& [QuantumFieldMax, PositionMax] = XMax;
      Field Tau = std::any_cast<Field>(PositionMax) - std::any_cast<Field>(PositionMin);
      int Sign = 1;
      while(Tau > 0){
	Tau -= beta;
	Sign *= -1;
      }
      while(Tau < -beta){
	Tau += beta;
	Sign *= -1;
      }
      return Sign*exp(Tau*mu)/(1+exp(-beta*mu));
    }
    
  };
    

}//namespace
