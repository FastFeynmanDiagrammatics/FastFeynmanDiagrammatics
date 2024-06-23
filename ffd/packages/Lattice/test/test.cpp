#include<ffd/std.hpp>
#include<ffd/core.hpp>
#include"../../PackagesMath.hpp"
#include"../../ClassTuple.hpp"
#include"../../PeriodicCoordinate.hpp"
#include"../../Lattice.hpp"

using namespace ffd::packages_math;
using namespace ffd::periodic_coordinate;
using namespace ffd::user_space;
using namespace ffd::lattice;


int main(){
  HoneycombLattice HoneyLattice({10, 10});
  for(auto j: {0,1}){
    for(auto k: {0,1}){
      std::cout<<HoneyLattice.BravaisLatticeObject.BravaisVectors[j][k]<<" ";
    }
    std::cout<<"\n";
  }
  auto X = CreateCoordinates(HoneyLattice);
  ffd::Get<2>(X)() = 1;
  ffd::Get<0>(X)() = 1;
  std::cout<<ffd::RealSpace(X)[0]<<" ddd "<<ffd::RealSpace(X)[1]<<std::endl;
  std::cout<<ffd::Get<0>(X)()<<" "<<ffd::Get<1>(X)()<<std::endl;
  for(auto yy: ffd::Get<0>(X).RealSpaceVectors){
    for(auto v: yy){
      std::cout<<v<<" ";
    }
    std::cout<<"\n";
  }
  for(auto yy: ffd::Get<1>(X).RealSpaceVectors){
    for(auto v: yy){
      std::cout<<v<<" ";
    }
    std::cout<<"\n";
  }

  std::cout<<"\n";
  for(int j=0; j<2; ++j){
    auto [low, upp] = HoneyLattice.BravaisLatticeObject.LowerUpperBounds[j];
    auto [x, y] = HoneyLattice.BravaisLatticeObject.BravaisVectors[j];
    std::cout<<"("<<low<<", "<<upp<<")\n";
    std::cout<<"bravais = ("<<x<<", "<<y<<")\n";
  }
}
