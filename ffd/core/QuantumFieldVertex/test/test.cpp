#include<ffd/std.hpp>
#include<ffd/core/QuantumField.hpp>
#include<ffd/core/QuantumFieldProduct_Dot.hpp>
#include<ffd/core/QuantumFieldVertex.hpp>


using namespace ffd::user_space;

struct Spacetime{
  double x;
  Spacetime(double x_): x(x_){}
};

struct Spacetime2{
  double x;
  Spacetime2(double x_): x(x_){}
};


int main(){
  Spacetime x1=1.23;
  Spacetime2 x2=2.34;
  auto GreenFunction = (Psi_(1))(x1)*(Bar(Psi_(1)))(x2);
  std::cout<<GreenFunction<<std::endl;
}
