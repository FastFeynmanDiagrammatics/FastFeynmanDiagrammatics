#include<ffd/std.hpp>
#include<ffd/core/CoreMath.hpp>
#include<ffd/core/QuantumField.hpp>
#include<ffd/core/QuantumFieldProduct_Dot.hpp>
#include<ffd/core/QuantumFieldVertex.hpp>
#include<ffd/core/QuantumFieldVertex.hpp>
#include<ffd/core/QuantumFieldGraph.hpp>
#include<ffd/core/FeynmanEdge.hpp>
#include<ffd/core/WickFunction.hpp>
#include<ffd/core/WickMatrix.hpp>


using namespace ffd::user_space;


int main(){
  double x = 1.;
  auto P_ = Psi_(1)*Psi_(-1);
  auto V_ = (Bar(P_)*P_)(x);
  auto G_ = V_|V_|V_;
  ffd::wick_function::WickFunction W_;
  auto Hubbard_minus = (Bar(Psi_(1))*Bar(Psi_(-1))*Psi_(1)*Psi_(-1))(x);
  W_ *= G_;

  ffd::wick_matrix::WickMatrix<double> DetM_(Psi_(1), 2);
  ffd::wick_matrix::WickMatrix<double> HafM_(Phi_(1), 3);
  ffd::wick_matrix::WickMatrix<double> PfaM_(Rho_(1), 3);
  for(int j=0; j < size(DetM_); ++j){
    for(int k=0; k < size(DetM_); ++k){
      DetM_(j, k, "assign") = std::sin(j)-cos(k*.2+3);
    }
  }
  for(int j=0; j < size(HafM_); ++j){
    for(int k=j; k < size(HafM_); ++k){
      HafM_(j, k, "assign") = std::sin(j)-cos(k*.2+3);
    }
  }

  for(int j=0; j < size(PfaM_); ++j){
    for(int k=j+1; k < size(PfaM_); ++k){
      PfaM_(j, k, "assign") = std::sin(j)-cos(k*.2+3);
    }
  }

  std::cout<<DetM_<<std::endl;
  std::cout<<HafM_<<std::endl;
  std::cout<<PfaM_<<std::endl;
  std::cout<<std::endl;

  std::cout<<DetM_(0, 0)*DetM_(1, 1)-DetM_(1, 0)*DetM_(0, 1)<<std::endl;
  std::cout<<ffd::wick_matrix::ComputePermutant(DetM_)<<std::endl;
  
}
