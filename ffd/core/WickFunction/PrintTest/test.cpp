#include<ffd/std.hpp>
#include<ffd/core/CoreMath.hpp>
#include<ffd/core/QuantumField.hpp>
#include<ffd/core/QuantumFieldProduct_Dot.hpp>
#include<ffd/core/QuantumFieldVertex.hpp>
#include<ffd/core/QuantumFieldVertex.hpp>
#include<ffd/core/QuantumFieldGraph.hpp>
#include<ffd/core/FeynmanEdge.hpp>
#include<ffd/core/WickFunction.hpp>


using namespace ffd::user_space;

std::optional<int> Find(const std::vector<int>& vec_, int val_){
  for(int j=0; j < std::size(vec_); ++j)if(vec_[j]==val_) return j; return {};
}


int main(){
  double x = 1.;
  auto P_ = Psi_(1)*Psi_(-1);
  auto V_ = (Bar(P_)*P_)(x);
  auto G_ = V_|V_|V_;
  ffd::wick_function::WickFunction W_;
  auto Hubbard_minus = (Bar(Psi_(1))*Bar(Psi_(-1))*Psi_(1)*Psi_(-1))(x);
  W_ *= G_;
  std::cout<<W_<<std::endl;
}





