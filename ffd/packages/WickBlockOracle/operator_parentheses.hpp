

namespace ffd::wick_block_oracle{

  char WickBlockOracle::operator()(QuantumField Q_) const{
    if(Q_.Dagger() < 0){
      Q_ = ffd::user_space::Bar(Q_);
    }
    auto which_block = PositionOfElement1InSetOfSets2(Q_, Blocks);
    assert( which_block.has_value() );
    return char( which_block.value() );
  }

}//namespace
