#include<ffd/std.hpp>
#include<ffd/core.hpp>
#include<ffd/tools.hpp>
#include<ffd/packages.hpp>


int main(){

  
  ffd::core_unit_test::UnitTest();
  std::cerr<<"ffd/core.hpp tests passed\n";

  
  auto time_test_tools = ffd::tools_unit_test::UnitTest();
  std::cerr<<"ffd/tools.hpp tests passed in "<<time_test_tools<<"s\n";

  
  auto time_test_packages = ffd::packages_unit_test::UnitTest();
  std::cerr<<"ffd/packages.hpp tests passed in "<<time_test_packages<<"s\n";

  
  return 0;
}
