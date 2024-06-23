#include<ffd/std.hpp>
#include<ffd/core.hpp>
#include"../../PackagesMath.hpp"
#include"../../ClassTuple.hpp"
#include"../../PeriodicCoordinate.hpp"

using namespace ffd::packages_math;
using namespace ffd::periodic_coordinate;
using namespace ffd::user_space;


int main(){
  PeriodicCoordinates<2, int, int, Real> X;
  ffd::Get<0>(X).RealSpaceVectors = {{1, 0}};
  ffd::Get<1>(X).RealSpaceVectors = {{0, 1}};
  ffd::Get<0>(X).Variable = -1;
  ffd::Get<1>(X).Variable = 12;
  ffd::Get<1>(X).LowerUpperBound = {0, 9};
  std::cout<<Sign(ffd::Get<1>(X))<<std::endl;
  auto Y = X;
  auto vec = ffd::RealSpace(Y);
  std::cout<<"print |vec| = "<<std::size(vec)<<std::endl;
  for(auto vv: vec){
    std::cout<<vv<<" ";
  }
  std::cout<<std::endl;
//   ffd::Get<1>(Y).RealSpaceMatrix = {{1,2}, {3.1,4}};
//   //  std::cout<<ffd::RealSpace(ffd::Get<0>(X))<<std::endl;
//   std::cout<<(ffd::Get<1>(X).RealSpaceMatrix.value()[1][0])<<std::endl;
}
