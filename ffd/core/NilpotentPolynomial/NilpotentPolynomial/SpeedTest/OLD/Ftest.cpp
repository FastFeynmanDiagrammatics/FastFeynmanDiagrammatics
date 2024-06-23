#include"../../../std.hpp"
#include"../../CoreMath.hpp"
#include"../../../tools/RandomDistributions.hpp"
#include"../../NilpotentPolynomial.hpp"
#include"../../../tools/Timer.hpp"


using namespace ffd::user_space;
using namespace std;
using ffd::random_distributions::RNGen;
using ffd::random_distributions::Proba;


constexpr int degree = 10;
constexpr long iterations = (1<<(23-degree));

template<int j>
using NilPoly = ffd::nilpotent_polynomial::NilpotentPolynomialF<j,Real>;


int main(){
  Real MaximalRelativeError1 = 0., MaximalRelativeError2 = 0.;
  Timer Timer_additions, Timer_multiplications, Timer_divisions, Timer_logarithm;
  Timer Timer_additions2, Timer_multiplications2, Timer_divisions2, Timer_logarithm2;
  cout<<degree<<"th degree polynomials, with random coefficients belonging to [-1,1]"<<endl;
   NilPoly<degree> X, Y, Z, W, K, J, R;
  //  NilPoly X{degree}, Y{degree}, Z{degree}, W{degree}, K{degree}, J{degree}, R{degree};
  // Z = "Zero";
  // Z *= 1.;

  // ffd::nilpotent_polynomial::ShiftByXi<10>(1.);

  
  //  cout<<setprecision(15)<<"The zero polynomial: "<<Z<<endl;
  Real maximal_error_ope = 0., relative_error_ope = 0.;
  for(long V=0; V< (1<<degree); ++V){
    X[V] = (2*Proba(RNGen)-1);
    Y[V] = (2*Proba(RNGen)-1);
  }

  Z = X*Y;

  // std::cerr<<X[V]<<std::endl;
  
  std::cerr<<"OKKKKKK";
  int xxx;


  // NilPoly<1> ZR = ffd::nilpotent_polynomial::ShiftByXi< (1<<degree) - 3 >(Z);

  // std::cout<<"ZR = "<<ZR<<std::endl;
  
  
  Timer_multiplications.ini();
  for(long iter=0; iter<iterations; ++iter){
    Z *= Y;
  }
  Timer_multiplications.fin();

  std::cout<<Z[(1<<degree)-1]<<std::endl;


  std::cout<<"Time multiplication = ";
  if(Timer_multiplications()/iterations > 1){
    std::cerr<<Timer_multiplications()/iterations<<" s"<<endl;
  }else if(Timer_multiplications()/iterations > 1.0e-3){
    std::cerr<<1'000*Timer_multiplications()/iterations<<" ms"<<endl;
  }else if(Timer_multiplications()/iterations > 1.0e-6){
    std::cerr<<1'000'000*Timer_multiplications()/iterations<<" \u03BCs"<<endl;
  }else{
    std::cerr<<1.0e9*Timer_multiplications()/iterations<<" ns"<<endl;
  }

  
  for(long V=0; V< (1<<degree); ++V){
    X[V] = (2*Proba(RNGen)-1);
    Y[V] = (2*Proba(RNGen)-1);
  }
  Z = X;

  Timer_divisions.ini();
  for(long iter=0; iter<iterations; ++iter){
    Z /= Y;
  }
  Timer_divisions.fin();
    std::cout<<Z[(1<<degree)-1]<<std::endl;
  //  std::cout<<Z<<std::endl;

  std::cout<<"Time division = ";
  if(Timer_divisions()/iterations > 1){
    std::cerr<<Timer_divisions()/iterations<<" s"<<endl;
  }else if(Timer_divisions()/iterations > 1.0e-3){
    std::cerr<<1'000*Timer_divisions()/iterations<<" ms"<<endl;
  }else if(Timer_divisions()/iterations > 1.0e-6){
    std::cerr<<1'000'000*Timer_divisions()/iterations<<" \u03BCs"<<endl;
  }else{
    std::cerr<<1.0e9*Timer_divisions()/iterations<<" ns"<<endl;
  }
  
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

  /*
  X[0] = Proba(RNGen);
  Y = ffd::nilpotent_polynomial::Logarithm(X);
  poly3 = poly1;
  for(long V=0; V < (1<<degree); ++V){
    poly3.coef[V] *= ffd::set_theory::cardinality_set(V, 1);
  }
  poly2 = poly3/poly1;
  for(long V=1; V < (1<<degree); ++V){
    poly2.coef[V] /= ffd::set_theory::cardinality_set(V);
  }
  poly2.coef[0] = log(poly1.coef[0]);
  poly5 = -poly1.coef[0];
    //  poly5.add(poly1, poly5);
    poly5 += poly1;
    poly5 *= 1/poly1.coef[0];
    poly6 = log(poly1.coef[0]);
    poly6 += poly5;
    //  cout<<poly1<<endl;
    //cout<<poly6<<endl;
    poly7 = poly5;
    for(int j=2; j<=degree; j++){
      poly8 = poly5*poly7;
      poly7 = poly8;
      poly8 *= ffd::math_space::MinusOneTo1(j+1)*1./j;
      poly6 += poly8;
      //    cout<<poly8<<endl;
      //cout<<poly6<<endl;
    }
    //  cout<<poly2<<endl;
    poly6 *= -1;
    poly6 = poly4 + poly6;
    poly2 *= -1;
    poly2 = poly4 + poly2;
    
  }
  
  cout<<"------LOGARITHM TEST-------"<<endl;
  cout<<"we have three ways of computing the logarithm:"<<endl;
  cout<<"||log^std P - log^2 P||_inf = "<<poly2.MaxAbsCoef()<<"    ||log^std P - log^2 P||/||log^std P|| = "<<poly2.MaxAbsCoef()/poly4.MaxAbsCoef()<<endl;
  //  cout<<"Error_{standard, 2}="<<poly2<<endl;
  cout<<"||log^std P - log^3 P||_inf = "<<poly6.MaxAbsCoef()<<"    ||log^std P - log^3 P||/||log^std P|| = "<<poly6.MaxAbsCoef()/poly4.MaxAbsCoef()<<endl;

  //  cout<<"Error_{standard, 3}="<<poly6<<endl;

  polya=poly10;
  Timer_additions.ini();
  for(int j=0; j<iterations; ++j){
    if(j%2==0)
      poly2 = poly10 + poly1;
    else
      poly10 = poly2 + poly1;
  }
  Timer_additions.fin();
  
  for(long V=0; V < (1<<degree); ++V){
    poly1.coef[V] = 2*(2*Proba(RNGen)-1);
    poly10.coef[V] = 2*(2*Proba(RNGen)-1);
  }

  
  for(long V=0; V < (1<<degree); ++V){
    poly1.coef[V] = 2*(2*Proba(RNGen)-1);
    poly10.coef[V] = 2*(2*Proba(RNGen)-1);
  }

  
  Timer_multiplications.ini();
  for(int j=0; j<iterations; ++j){
    if(j%2==0)
      poly2 = poly10*poly1;
    else
      poly10 = poly2*poly1;
  }
  Timer_multiplications.fin();


  for(long V=0; V < (1<<degree); ++V){
    poly1.coef[V] = 2*(2*Proba(RNGen)-1);
    poly10.coef[V] = 2*(2*Proba(RNGen)-1);
  }

  
  Timer_divisions.ini();
  for(int j=0; j<iterations; ++j){
    if(j%2==0)
      poly2 = poly10/poly1;
    else
      poly10 = poly2/poly1;
  }
  Timer_divisions.fin();

  for(long V=0; V < (1<<degree); ++V){
    poly1.coef[V] = 2*(2*Proba(RNGen)-1);
    poly10.coef[V] = 2*(2*Proba(RNGen)-1);
  }


  Timer_logarithm.ini();
  poly1.coef[0] = 2.231323;
  poly3 = poly1;
  for(int j=0; j<iterations; ++j){
    if(j%2==0)
      poly2.logarithm(poly3);
    else
      poly3.logarithm(poly2);
  }
  Timer_logarithm.fin();
    */  
}
