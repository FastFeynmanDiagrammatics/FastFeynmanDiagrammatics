#include"../../../std.hpp"
#include"../../CoreMath.hpp"
#include"../../../tools/RandomDistributions.hpp"
#include"../../NilpotentPolynomial.hpp"
#include"../../../tools/Timer.hpp"


using namespace ffd::user_space;
using namespace std;
using ffd::random_distributions::RNGen;
using ffd::random_distributions::Proba;


const int nilpot = 1, degree = 10;
const long iterations = 1<<(24-degree);

template<int j>
using NilPoly = ffd::nilpotent_polynomial::NilpotentPolynomialF<j>;


int main(){
  Real MaximalRelativeError1 = 0., MaximalRelativeError2 = 0.;
  Timer Timer_additions, Timer_multiplications, Timer_divisions, Timer_logarithm;
  Timer Timer_additions2, Timer_multiplications2, Timer_divisions2, Timer_logarithm2;
  cout<<degree<<"th degree polynomials, with random coefficients belonging to [-1,1]"<<endl;
  NilPoly<degree> X, Y, Z;//, W, K, J, R;
  Z = "Zero";
  Z *= 1.;

  for(ffd::BinaryInt V=0; V<(1<<3); ++V){
    for(ffd::BinaryInt S=V; S!=0; S=((S-1)&V)){
      std::cout<<S<<" "<<V<<" "<<std::endl;
    }
  }
  
  Real maximal_error_ope = 0., relative_error_ope = 0.;
  for(long V=0; V< (1<<degree); ++V){
      X[V] = (2*Proba(RNGen)-1);
      Y[V] = 2.7*(2*Proba(RNGen)-1);
  }
  Z = X;

  Timer_multiplications.ini();
  for(long iter=0; iter<iterations; ++iter){
    Z *= Y;
  }
  Timer_multiplications.fin();
  std::cout<<Z[(1<<degree)-1]<<std::endl;

  std::cerr<<"CPU time multiplication = ";
  if(Timer_multiplications()/iterations > 1){
    std::cerr<<Timer_multiplications()/iterations<<" s"<<endl;
  }else if(Timer_multiplications()/iterations > 1.0e-3){
    std::cerr<<1'000*Timer_multiplications()/iterations<<" ms"<<endl;
  }else if(Timer_multiplications()/iterations > 1.0e-6){
    std::cerr<<1'000'000*Timer_multiplications()/iterations<<" \u03BCs"<<endl;
  }else{
    std::cerr<<1.0e9*Timer_multiplications()/iterations<<" ns"<<endl;
  }

  Z = X;
  Timer_divisions.ini();
  for(long iter=0; iter<iterations; ++iter){
    Z /= Y;
  }
  Timer_divisions.fin();
  //  std::cout<<Z<<std::endl;
  std::cout<<Z[(1<<degree)-1]<<std::endl;
  
  std::cerr<<"CPU time division       = ";
  if(Timer_divisions()/iterations > 1){
    std::cerr<<Timer_divisions()/iterations<<" s"<<endl;
  }else if(Timer_divisions()/iterations > 1.0e-3){
    std::cerr<<1'000*Timer_divisions()/iterations<<" ms"<<endl;
  }else if(Timer_divisions()/iterations > 1.0e-6){
    std::cerr<<1'000'000*Timer_divisions()/iterations<<" \u03BCs"<<endl;
  }else{
    std::cerr<<1.0e9*Timer_divisions()/iterations<<" ns"<<endl;
  }
  std::cout<<Z[(1<<degree)-1]<<std::endl;
  /*

  
  for(long iter=0; iter<iterations; ++iter){
    for(long V=0; V< (1<<degree); ++V){
      X[V] = (2*Proba(RNGen)-1);
      Y[V] = (2*Proba(RNGen)-1);
      Z[V] = (2*Proba(RNGen)-1);
      W[V] = (2*Proba(RNGen)-1);
    }
    K = (X+Z)/(Y-W);
    J = (Y-W)/(X+Z);
    R = K*J + J*K;
    R /= 2.;
    R -= 1.;
    if(maximal_error_ope < R.MaxAbsCoef()){
      maximal_error_ope = R.MaxAbsCoef();
    }
    Real RelErr = 2.*R.MaxAbsCoef()/(J.MaxAbsCoef()+K.MaxAbsCoef());
    if(relative_error_ope < RelErr){
      relative_error_ope = RelErr;
    }
  }
  cout<<"------RING OPERATIONS TEST-------"<<endl;
  cout<<"this should be zero (mathematically):"<<endl;
  cout<<"Maximal ||((X+Z)/(Y-W))*((Y-W)/(X+Z)) - 1||_inf = "<<maximal_error_ope<<endl;
  cout<<"Maximal Relative Error = "<<relative_error_ope<<endl;
  */
}
