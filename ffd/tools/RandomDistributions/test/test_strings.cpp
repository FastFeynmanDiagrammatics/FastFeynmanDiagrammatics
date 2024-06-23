#include<ffd/std.hpp>
#include<ffd/core.hpp>
#include<ffd/tools.hpp>
#include<ffd/packages.hpp>

int main(){

  for (uint i=0; i<100; ++i){
    std::cout << ffd::random_distributions::RandomString(10) << std::endl;
  }

}
