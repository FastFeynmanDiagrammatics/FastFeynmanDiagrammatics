

namespace ffd::richardson_extrapolation{

  void rich_extrap_three_points(){
    using ffd::core_math::Pi;
    using namespace std;

    const int order = 1;
    std::vector<Complex> sequence(1+2*order);
    const Real Beta = 2, tau = 0.1, xi = .2;
    auto omega_n = [Beta](int n)->Real{return 2.*Pi*(n+.5)/Beta;};
    auto summand = [omega_n, tau, xi, Beta](int n)
      ->Complex{return (1/Beta)*exp(Complex(0, -omega_n(n)*tau))*
	(1./(Complex(0, omega_n(n))-xi)-1./Complex(0, omega_n(n)));};
    sequence[0] = summand(0)+summand(-1);
    for(int u=1; u < 1+2*order; ++u){
      sequence[u] = sequence[u-1];
      for(int j = 1<<(u-1); j < 1<<u; ++j){
	sequence[u] += summand(j) + summand(-j-1); 
      }
    }
    //    std::cerr<<summand(0)<<" "<<sequence[0]<<" "<<sequence[1]<<" "<<sequence[2]<<" "<<sequence[3]<<
    //      " "<<sequence[4]<<std::endl;
    auto [value, error] = RichardsonExtrapolationThreePoints(sequence);
    std::cout<<sequence[2*order]<<" "<<std::setprecision(20)<<value<<" "<<error<<std::endl;
    std::cout<<(-exp(-xi*tau)/(1+exp(-Beta*xi))+.5)<<std::endl;
  }

}//namesapce
