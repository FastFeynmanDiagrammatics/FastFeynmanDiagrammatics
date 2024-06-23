

namespace ffd::user_space::imaginary_time_lattice_propagator{

  Real

  Constrain1BetweenZeroAnd2ReturnSign(Real& tau,
				      Real beta,
				      bool is_fermion){
    Real sign = 1;

    while( tau < 0 ){
      tau += beta;
      if(is_fermion){
	sign *= -1;
      }
    }

    while( tau > beta){
      tau -= beta;
      if(is_fermion){
	sign *= -1;
      }
    }


    return sign;
  }

}//namespace
