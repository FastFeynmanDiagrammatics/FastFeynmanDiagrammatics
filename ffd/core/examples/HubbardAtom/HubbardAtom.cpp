#include<ffd/core.hpp>
#include<ffd/tools/RandomDistributions.hpp>
#include<ffd/tools/Determinant.hpp>
#include<ffd/packages/Permutant.hpp>
#include<ffd/tools/AnalyticalDerivative.hpp>
#include<ffd/tools/Timer.hpp>
#include"BarePropagator.hpp"
#include"AnalyticalFunctions.hpp"


//we use the namespace containing declarations commonly used by a ffd user
using namespace ffd::user_space;
using namespace ffd::user_space::core_examples;


//inverse temperature
const Real Beta = 1;
//Chemical potential for U=0
const Real mu0 = -1;
constexpr bool PrintToScreen = true;


//the first argument is the order, the second is the number of iterations
int main(int argc, char* argv[]){
  if(argc<2){std::cerr<<"!!!order and iterations missing\n"; return 1;}
  if(argc<3){std::cerr<<"!!!iterations missing\n"; return 1;}
  const int order = std::stoi(argv[1]);
  const long iterations = std::stol(argv[2]);
  std::cout<<"order = "<<order<<", MC_steps = "<<iterations<<std::endl;

  //we start defining the interaction part of the hamiltonian
  //and the external vertices
  //Psi_(s) is a  QuantumFieldProduct composed of one
  //fermionic QuantumField of component s
  //Rho_(s) would be a Majorana
  //Eta_(s) would be a boson
  //Phi_(s) would be a scalar particle
  std::map<int, decltype(Psi_(1))> n;
  for(auto s: {-1, 1}){
    auto Psi = Psi_(s);
    //n[s] is our density QuantumFieldProduct
    n[s] = Bar(Psi)*Psi;
  }
  //we define the imaginary time in the most elementary way possible
  //we leave more advanced definitions outside core
  Real tau = .5*Beta;
  //HubbardVertex is a QuantumFieldDot (that can be converted to QuantumFieldVertex)
  auto HubbardVertex = (n[1]*n[-1])(tau);
  QuantumFieldGraph Graph;
  for(int j=0; j < order; ++j){
    //the bitwise or operator | is used to add a QuantumFieldVertex to a QuantumFieldGraph
    Graph |= HubbardVertex;
  }
  Real tau1 = 0.5, tau2 = 0.5;
  //we define a GreenFunction QuantumFieldVertex
  auto GreenFunctionVertex = Bar(Psi_(1))(tau1)*Psi_(1)(tau2);
  //the Graph contains all possible external vertices
  //here we have a Green Function external vertex and
  //an interaction potential external vertex
  Graph |= GreenFunctionVertex;
  Graph |= HubbardVertex;
  std::cout<<Graph<<std::endl;

  //we introduce the polynomials in which we will store the
  //sum of Feynman diagrams
  NilPolynomial Z(order), G(order), D(order);
  //these are accumulation variables for the Monte Carlo
  Real G_acc = 0, Z_acc = 0, G2_acc = 0, Z2_acc = 0, D_acc = 0, D2_acc = 0;
  Timer MapClock, SubsetsClock, SubsetsZClock, SubsetsGClock, RecursionClock, DetClock, MatrixClock, SignZClock, PreWickClock;
  PreWickClock.ini();
  std::vector<short> SignWZ(1<<order);
  std::vector<std::vector<ffd::wick_matrix::WickMatrix<FeynmanEdge>>> EdgeMatrices(1<<order);
  for(BinaryInt S=0; S < 1<<order; ++S){
    //      WickClock.ini();
    //we fill the WickFunction with QuantumFieldVertex elements
    WickFunction WZ;
    for(auto j: VectorOfBinaryDigitsOf(S)){
      WZ *= Graph(j);
    }
    SignWZ[S] = WZ.Sign;
    //      WickClock.fin();
    EdgeMatrices[S] = CreateFeynmanEdgeWickMatrices(WZ);      
  }
  PreWickClock.fin();

  //MC loop
  for(long t_MC=0; t_MC < iterations; ++t_MC){
    //random initialization (direct sampling MC)
    for(int j=0; j < order; ++j){
      Graph[j][0].Position = Beta*Proba(RNGen);
    }

    //a FeynmanEdgeMap associate to each FeynmanEdge in the Graph
    //a value (the type of this value is chosen to be Real in this
    //case, but one can use for example polynomials)
    //We give an std::function representing the bare propagator as argument
    MapClock.ini();
    auto EdgeMap = CreateFeynmanEdgeMap<Real>(Graph, BuildBarePropagator(Beta, mu0));
    MapClock.fin();
    
    // if(t_MC == 0 && PrintToScreen){
    //   std::cout<<EdgeMap<<std::endl;
    // }

    
    SubsetsClock.ini();
    SubsetsZClock.ini();
    //looping over subsets of 11...111
    for(BinaryInt S=0; S < (1<<order); ++S){
      //computing the partition function
      //The WickFunction sorts and stores the QuantumFields
      //The fields are stored in different correlation blocks,
      //which can be imposed by giving an appropriate std::function
      //when constructing the object. Default blocks are
      //by component and statistics
      //      WickClock.ini();
      //      WickClock.fin();
      if(t_MC == 0 && PrintToScreen){
	std::cout<<"S = "<<std::bitset<sizeof(S)*4>(S)<<std::endl;
      }
      //once the WickFunction is filled, we can compute
      //a WickMatrix for each correlation block
      MatrixClock.ini();
      auto const& MatricesS = EdgeMatrices[S];
      const int SizeEdgeMatricesS = std::size(MatricesS);
      std::vector<ffd::wick_matrix::WickMatrix<Real>> Matrices(SizeEdgeMatricesS);
      for(int j=0; j < SizeEdgeMatricesS; ++j){
	Matrices[j] = CreateNumericWickMatrix(MatricesS[j], EdgeMap);
      }
      MatrixClock.fin();
      //the weight of the vertices is (-1)^{number of vertices} for the
      //Hubbard model. We also multiply by the Sign of the WickFunction
      //which takes into account the anticommutation of operators
      //      SignZClock.ini();
      Z[S] = std::pow(-1, ffd::set_theory::CardinalitySet(S))*SignWZ[S];
      //      SignZClock.fin();
      //we can then compute the determinants of the WickMatrix
      //In order to be general, we use a ``permutant'' function
      //that will be able to handle determinants, pfaffians, etc..
      //the user must furnish at least one function ffd::determinant::Determinant
      //(we use the one from ffd/tools/Determinant.hpp)
      DetClock.ini();
      for(auto& M: Matrices){
	Z[S] *= ComputePermutant(M);
      }
      DetClock.fin();
    }
    SubsetsZClock.fin();
    
    SubsetsGClock.ini();
    for(BinaryInt S=0; S < (1<<order); ++S){    
      //computing the (disconnected) Green's function (the density in this case)
      WickFunction WA;
      for(auto j: VectorOfBinaryDigitsOf(S)){
	WA *= Graph(j);
      }
      WA *= Graph(order);
      if(t_MC == 0 && PrintToScreen){
	std::cerr<<WA<<std::endl;
      }
      auto EdgeMatrices = CreateFeynmanEdgeWickMatrices(WA);
      std::vector<ffd::wick_matrix::WickMatrix<Real>> Matrices;
      for(auto const& EdgeMatrix: EdgeMatrices){
	Matrices.push_back(CreateNumericWickMatrix(EdgeMatrix, EdgeMap));
      }

      G[S] = std::pow(-1, ffd::set_theory::CardinalitySet(S))*WA.Sign;
      for(auto& M: Matrices){
	if(t_MC == 0 && PrintToScreen){
	  std::cout<<M<<std::endl;
	}
	G[S] *= ComputePermutant(M);
      }
    }
    SubsetsGClock.fin();

    
    for(BinaryInt S=0; S < (1<<order); ++S){
      //computing the (disconnected) interaction energy
      WickFunction WAA;
      for(auto j: VectorOfBinaryDigitsOf(S)){
	WAA *= Graph(j);
      }
      WAA *= Graph(order+1);
      if(t_MC == 0 && PrintToScreen){
	std::cerr<<WAA<<std::endl;
      }
      auto EdgeMatrices = CreateFeynmanEdgeWickMatrices(WAA);
      std::vector<ffd::wick_matrix::WickMatrix<Real>> Matrices;
      for(auto const& EdgeMatrix: EdgeMatrices){
	Matrices.push_back(CreateNumericWickMatrix(EdgeMatrix, EdgeMap));
      }
      
      D[S] = std::pow(-1, ffd::set_theory::CardinalitySet(S))*WAA.Sign;
      for(auto& M: Matrices){
	if(t_MC == 0 && PrintToScreen){
	  std::cout<<M<<std::endl;
	}
	D[S] *= ComputePermutant(M);
      }
      
    }//end of loop over subsets
    SubsetsClock.fin();

    RecursionClock.ini();
    //eliminating disconnected diagrams
    G /= Z;
    RecursionClock.fin();
    D /= Z;
    
    //accumulate G, D and Z
    Real G_temp = G[(1<<order)-1]/ffd::core_math::Factorial(order);
    Real D_temp = D[(1<<order)-1]/ffd::core_math::Factorial(order);
    Real Z_temp = Z[(1<<order)-1]/ffd::core_math::Factorial(order);
    G_acc += G_temp;
    G2_acc += std::pow(G_temp, 2);
    D_acc += D_temp;
    D2_acc += std::pow(D_temp, 2);
    Z_acc += Z_temp;
    Z2_acc += std::pow(Z_temp, 2);
  }
  G_acc /= iterations; G2_acc /= iterations;   D_acc /= iterations; D2_acc /= iterations; Z_acc/= iterations; Z2_acc /= iterations;
  G2_acc -= std::pow(G_acc, 2);   D2_acc -= std::pow(D_acc, 2); Z2_acc -= std::pow(Z_acc, 2);
  G2_acc = std::sqrt(G2_acc/iterations);  D2_acc = std::sqrt(D2_acc/iterations); Z2_acc = std::sqrt(Z2_acc/iterations);
  Real RadiusAnalytical = 1;
  std::cout<<"n^{"<<order<<"} = "<<G_acc<<" +/- "<<G2_acc<<" (Monte Carlo)"<<std::endl;
  ffd::analytical_derivative::AnalyticalDerivative n_analytical(BuildDensity(Beta, mu0));
  std::cout<<"n^{"<<order<<"} = "<<std::real(n_analytical.CoefficientOfzTo1(order, 0., RadiusAnalytical))<<" (Exact)"<<std::endl;
  std::cout<<"D^{"<<order<<"} = "<<D_acc<<" +/- "<<D2_acc<<" (Monte Carlo)"<<std::endl;
  ffd::analytical_derivative::AnalyticalDerivative D_analytical(BuildDoubleOccupancy(Beta, mu0));
  std::cout<<"D^{"<<order<<"} = "<<std::real(D_analytical.CoefficientOfzTo1(order, 0., RadiusAnalytical))<<" (Exact)"<<std::endl;
  std::cout<<"Z^{"<<order<<"} = "<<Z_acc<<" +/- "<<Z2_acc<<" (Monte Carlo)"<<std::endl;
  ffd::analytical_derivative::AnalyticalDerivative Z_analytical(BuildPartitionFunction(Beta, mu0));
  std::cout<<"Z^{"<<order<<"} = "<<std::real(Z_analytical.CoefficientOfzTo1(order, 0., RadiusAnalytical)/Z_analytical.CoefficientOfzTo1(0,0., RadiusAnalytical))<<" (Exact)"<<std::endl;
  std::cout<<"Time Pre-Wick = "<<1e3*PreWickClock()<<" m s"<<std::endl;
  std::cout<<"Time Map = "<<1e3*MapClock()/iterations<<" m s"<<std::endl;
  std::cout<<"Time Subsets/steps = "<<1e3*SubsetsClock()/iterations<<" m s"<<std::endl;
  std::cout<<"Time Subsets Z /steps = "<<1e3*SubsetsZClock()/iterations<<" m s"<<std::endl;
  std::cout<<"Time Sign Z/steps = "<<1e3*SignZClock()/iterations<<" m s"<<std::endl;
  std::cout<<"Time WickMatrix Z/steps = "<<1e3*MatrixClock()/iterations<<" m s"<<std::endl;
  std::cout<<"Time Determinant/steps = "<<1e3*DetClock()/iterations<<" m s"<<std::endl;
  std::cout<<"Time Subsets G /steps = "<<1e3*SubsetsGClock()/iterations<<" m s"<<std::endl;
  std::cout<<"Time Recursion/steps = "<<1e3*RecursionClock()/iterations<<" m s"<<std::endl;
}

