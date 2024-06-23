

namespace ffd::wick_block_oracle{

  template<typename Field>
  WickBlockOracle::WickBlockOracle(ffd::user_space::QuadraticAction<Field> const& A_){
    for(auto const& term: A_){
      AddQuadraticActionTerm(term);
    }
  }

}//namespace
