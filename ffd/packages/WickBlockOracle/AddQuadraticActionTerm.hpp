

namespace ffd::wick_block_oracle{

  template<typename Field>
  void WickBlockOracle::AddQuadraticActionTerm(ffd::user_space::QuadraticActionTerm<Field>
					       const& term_){
    assert( std::size(term_) == 2);
    std::array<QuantumField, 2> correlated_fields;
    for(auto j: {0, 1}){
      auto const& quantum_field_dot = term_[j];
      auto first_and_only_field = quantum_field_dot[0];
      if( first_and_only_field.Dagger() < 0){
	first_and_only_field = ffd::user_space::Bar(first_and_only_field);
      }
      correlated_fields[j] = first_and_only_field;
    }
    AddRelation1ToRelationSets2(correlated_fields, Blocks);
  }


}//namespace

