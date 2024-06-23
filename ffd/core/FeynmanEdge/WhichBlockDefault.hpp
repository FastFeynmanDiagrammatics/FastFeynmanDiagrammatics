

namespace ffd::feynman_edge{

  char WhichBlockDefault(ffd::quantum_field::QuantumField Q_){
    return Q_.Component() + 31*( 1 - Q_.Dagger()*Q_.Dagger() + 2*( !Q_.IsFermion() ) );
  }

}//namespace 
