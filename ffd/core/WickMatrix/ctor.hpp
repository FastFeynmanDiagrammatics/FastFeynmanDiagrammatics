namespace ffd::wick_matrix{

  template<typename Field, typename intType>
  WickMatrix<Field, intType>::WickMatrix(QuantumField const& Q_,
					 intType NumberDestructionOperators_):
    NumberDestructionOperators(NumberDestructionOperators_),
    IsFermion( Q_.IsFermion() ), NotHermitian( Q_.Dagger() != 0 ){
    intType const& n_ = NumberDestructionOperators;
    if(NotHermitian){
      Components.resize(n_*n_);
    }else{
      if(IsFermion){
	Components.resize((n_*(n_-1))/2);
      }else{
	Components.resize((n_*(n_+1))/2);
      }
    }
  }

  
  template<typename Field, typename intType>
  WickMatrix<Field, intType>::WickMatrix(const QuantumFieldProduct& Q_,
				intType NumberDestructionOperators_):
    WickMatrix(Q_[0], NumberDestructionOperators_) {}


}//namespace ffd::wick_matrix
