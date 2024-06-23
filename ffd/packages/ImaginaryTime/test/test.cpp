#include<ffd/std.hpp>
#include<ffd/core.hpp>
#include"../../PackagesMath.hpp"
#include"../../ClassTuple.hpp"
#include"../../PeriodicCoordinate.hpp"
#include"../../ImaginaryTime.hpp"

using namespace ffd::user_space;

int main(){
  ffd::imaginary_time::ImaginaryTime ImagTime{2.2, 0};
  auto tau = CreatePeriodicCoordinates(ImagTime);
  std::cout<<ffd::Get<0>(tau).Variable<<std::endl;
}
