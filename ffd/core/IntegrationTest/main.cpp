#include<ffd/core.hpp>
#include"BarePropagator.hpp"
#include"ShufflePermutation.hpp"
#include"MonteCarloDeterminant.hpp"
#include"../../tools/Determinant.hpp"

using namespace ffd::user_space;
using namespace ffd::core_math;
using namespace ffd::core_integration_test;
using namespace ffd::gauss_determinant;


using LReal = long double;


const LReal Beta = 1.l;
const LReal Mu0 = -0.25l;
// const long samples = 1<<25;
// const long shuffles = 0;


int main(int argc, char* argv[]){
  if(argc < 2){
    throw std::invalid_argument("give max_order argument in command line");
  }
  const int max_order = std::stoi(argv[1]);
  std::map<char, decltype(Psi_(1))> n;
  for(char s: {-1, 1}){
    n[s] = Bar(Psi_(s))*Psi_(s);
  }
  LReal tau = .5;
  auto HubbardVertex = (n[1]*n[-1])(tau);
  for(int order=2; order <= max_order; ++order){

    std::uniform_real_distribution<long double> Proba(0.l, 1.l);
    QuantumFieldGraph G = n[1](tau);
    for(int j=0; j<order; ++j){
      HubbardVertex.Position = Beta*Proba(gen);//Factorial(j)*132.12312*Beta;
      G |= HubbardVertex;
    }

    WickFunction W;
    for(int j=0; j<=order; ++j){
      W *= G(j);
    }

    auto EdgeMatrices = CreateFeynmanEdgeWickMatrices(W);

    auto G0 = BarePropagator(Beta, Mu0);
    
    auto EdgeMap = CreateFeynmanEdgeMap<LReal>(G, G0);  

    LReal ffd_result = 1.l;
    for(auto const& EdgeMatrix: EdgeMatrices){
      auto NumericMatrix = CreateNumericWickMatrix(EdgeMatrix, EdgeMap);
      ffd_result *= Determinant(NumericMatrix);
    }
    ffd_result /= Factorial(order)*std::pow(-1, order);

    LReal exact_result = 0, Z0 = 0, alpha = 1/(1+exp(-Beta*Mu0));
    for(auto n_up: {0, 1}){
      for(auto n_down: {0, 1}){
	exact_result += (n_up-alpha)*exp(Beta*Mu0*(n_up+n_down))*
	  std::pow(-Beta*(n_up - alpha)*(n_down - alpha), order)/Factorial(order);
	Z0 += exp(Beta*Mu0*(n_up+n_down));
      }
    }
    exact_result /= Z0;
		    
    assert(std::abs(ffd_result-exact_result)/std::abs(exact_result) < 1000*std::numeric_limits<LReal>::epsilon());
  }
}
