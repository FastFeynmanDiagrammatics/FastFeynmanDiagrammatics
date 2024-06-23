#include<ffd/std.hpp>
#include<ffd/core.hpp>
#include"../../PackagesMath.hpp"
#include"../../QuadraticAction.hpp"
#include"../../WickBlockOracle.hpp"

using namespace ffd::user_space;

int main(){
  std::array<int, 2> x1{0,1}, x2{-1, 0};
  QuadraticAction<Real> H_0 = (Bar(Psi_(1))(x1)*Psi_(1)(x2))*(-1.);
  H_0 += FlipSpin(H_0);
  auto Oracle = ffd::wick_block_oracle::CreateWickBlockOracleFromAction1(H_0);
  std::cout<<Oracle(Bar(Psi_(1))[0])<<" "<<Oracle(Psi_(-1)[0])<<" "<<std::endl;
}
