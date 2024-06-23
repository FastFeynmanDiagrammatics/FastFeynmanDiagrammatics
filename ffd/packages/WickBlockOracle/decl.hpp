namespace ffd::wick_block_oracle{

  using namespace ffd::relation_sets;
  using ffd::quantum_field::QuantumField;

  
  class WickBlockOracle{
  public:
    RelationSets<QuantumField> Blocks;

    
    template<typename Field>
    WickBlockOracle(ffd::user_space::QuadraticAction<Field> const&
		    quadratic_action_encoding_the_non_interacting_correlations);

    
    template<typename Field>
    void AddQuadraticActionTerm(ffd::user_space::QuadraticActionTerm<Field> const&
				a_quadratic_term_of_the_quadratic_action);

    
    char
    operator()(
	       QuantumField
	       a_quantum_field_to_which_we_associate_a_correlation_block_stored_as_a_numeric_char
	       )
      const;
    
  };

  
  
  auto size(WickBlockOracle const& W_){
    return size(W_.Blocks);
  }

  
}//namespace
