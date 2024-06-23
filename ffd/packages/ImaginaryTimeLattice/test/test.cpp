#include<ffd/std.hpp>
#include<ffd/core.hpp>
#include"../../PackagesMath.hpp"
#include"../../ClassTuple.hpp"
#include"../../PeriodicCoordinate.hpp"
#include"../../Lattice.hpp"
#include"../../ImaginaryTime.hpp"
#include"../../ImaginaryTimeLattice.hpp"

using namespace ffd::user_space;

const Real Beta = 2.2;
const int L = 10;

int main(){
  ffd::imaginary_time::ImaginaryTime ImagTime{Beta};
  auto tau = CreateCoordinates<2>(ImagTime);
  ffd::lattice::HoneycombLattice H({L, L});
  auto x = CreateCoordinates(H);
  auto y = x;
  auto Y = tau*x;
  auto X = Y;
  ffd::Get<1>(X)() = 6;
  std::cout<<ffd::Get<0>(X)()<<std::endl;
  std::cout<<ffd::Get<1>(X)()<<std::endl;
  std::cout<<ffd::Get<2>(X)()<<std::endl;
  std::cout<<ffd::Get<3>(X)()<<std::endl;
  std::cout<<ffd::RealSpace(X)[0]<<std::endl;
  std::cout<<ffd::RealSpace(X)[1]<<std::endl;
}
