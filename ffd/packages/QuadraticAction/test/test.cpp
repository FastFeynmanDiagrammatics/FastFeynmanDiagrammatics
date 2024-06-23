#include<ffd/std.hpp>
#include<ffd/core.hpp>
#include"../../QuadraticAction.hpp"

using namespace ffd::user_space;

int main(){
  Real x1 = 0, x2= 1;
  RealQuadraticAction H_0 = (Psi_(-1)(x1)*Psi_(1)(x2))*1.;
  H_0 += Bar(H_0);
  H_0 += FlipSpin(H_0);
  std::cout<< H_0 << "\n";
  std::cout<<sizeof(float)<<std::endl;
  
}
