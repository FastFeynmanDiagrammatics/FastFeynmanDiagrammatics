

namespace ffd::user_space::imaginary_time_lattice_propagator{

  Real
  GreenFunctionZetaEnergyTau(bool is_fermion,
			     Real beta,
			     Real E,
			     Real tau){
    return is_fermion ?
      exp( (.5*beta - tau)*E )/cosh( .5*beta*E ):
      -exp( (.5*beta - tau)*E )/sinh( .5*beta*E);
  }

}//namespace
