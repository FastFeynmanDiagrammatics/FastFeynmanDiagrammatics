

namespace ffd::richardson_extrapolation{

  void rich_extrap_integer(){
    using ffd::core_math::Pi;
    using namespace std;

    const int order = 10;
    std::vector<Complex> sequence(order+1);
    const Real Beta = 100, tau = .1213, xi = .3213;
    auto omega_n = [Beta](int n)->Real{return 2.*Pi*(n+.5)/Beta;};
    auto summand = [omega_n, tau, xi, Beta](int n)
      ->Complex{return (1/Beta)*exp(Complex(0, -omega_n(n)*tau))*
	(1./(Complex(0, omega_n(n))-xi)-1./Complex(0, omega_n(n)));};
    sequence[0] = summand(0);
    for(int u=1; u <= order; ++u){
      sequence[u] = sequence[u-1];
      for(int j = 1<<(u-1); j < 1<<u; ++j){
	sequence[u] += summand(j);
      }
    }
    //    std::cerr<<summand(0)<<" "<<sequence[0]<<" "<<sequence[1]<<" "<<sequence[2]<<" "<<sequence[3]<<
    //      " "<<sequence[4]<<std::endl;
    Real phase = 0*(-2*Pi*tau/Beta);
    auto values = RichardsonExtrapolationIntegerPowersPhase(sequence, phase);
    Real exact = (-exp(-xi*tau)/(1+exp(-Beta*xi))+.5);
    std::cout<<std::setprecision(20)<<0<<" "<<sequence[order]*2.-exact<<std::endl;
    for(std::size_t j=0; j < size(values); ++j){
      std::cout<<j+1<<" "<<std::setprecision(20)<<values[j]*2.-exact<<std::endl;
    }
    std::cout<<exact<<std::endl;

    // Complex summ = -exact/2+summand(0);
    // const int order_R = 10;
    // std::vector<Complex> summ_R(order_R, summ);
    // Real epsilon_R = .9;
    // int n_max = 998;
    // for(int j=1; j < n_max; ++j){
    //   summ += summand(j);
    //   for(int k=0; k<order_R; ++k){
    // 	summ_R[k] += summand(j)*pow(1.-epsilon_R/(1<<k), j);
    //   }
    //   std::cout<<j<<" "<<std::real(summ)<<" "<<std::real(summ_R[0])<<" "<<std::real(summand(j))*j*j<<std::endl;
    // }
    // for(int k=0; k<std::size(summ_R); ++k){
    //   std::cout<<k<<" "<<summ_R[k]<<std::endl;
    // }
    // auto values_R = RichardsonExtrapolationIntegerPowers(summ_R);
    // for(int k=0; k<std::size(values_R); ++k){
    //   std::cout<<k<<" "<<values_R[k]<<std::endl;
    // }

    
    // Complex summ = -exact/2+summand(0);
    // const int order_R = 20;
    // std::vector<Complex> summ_R(order_R, summ);
    // std::vector<Real> radius_vec(order_R);
    // Real a = 0, b= 1.0001;
    // for(int j=0; j<order_R; ++j){
    //   radius_vec[j] = .5*(a+b)+(b-a)*cos(Pi*(j+.5)/order_R)/2;
    // }
    // int n_max = 100;
    // for(int j=1; j < n_max; ++j){
    //   summ += summand(j);
    //   for(int k=0; k<order_R; ++k){
    // 	summ_R[k] += summand(j)*pow(radius_vec[k], j);
    //   }
    //   std::cout<<j<<" "<<std::real(summ)<<" "<<std::real(summ_R[0])<<" "<<std::real(summand(j))*j*j<<std::endl;
    // }
    // std::vector<Real> summ_R_real;
    // for(int k=0; k<std::size(summ_R); ++k){
    //   std::cout<<k<<" "<<summ_R[k]<<std::endl;
    //   //      summ_R_real.push_back(summ_R[k]);
    // }
    // ffd::chebyshev_polynomial::ChebyshevPolynomial<Complex> Poly(summ_R, {0, 1});
    // std::cout<<Poly(1)<<" "<<Poly.ErrorEstimate()<<std::endl;

    //    Real Radius_R = 1;
    auto summand2 = [](int n)->Real{return 1/pow(n+1., 2);};
    exact = 2*Pi*Pi/6;
    
    Complex summ = 0.;
    const int order_cheby = 2000;
    std::vector<Complex> cheby_values(2*order_cheby);
    for(int j=1; j <= order_cheby; ++j){
      summ += summand2(j);
      cheby_values[order_cheby-j] = summ;
      //      std::cout<<order_cheby-j<<" "<<summ<<std::endl;
    }
    for(int j=0; j<order_cheby; ++j){
      cheby_values[order_cheby+j] = cheby_values[order_cheby-j-1]*(-1.);
    }
    ffd::chebyshev_polynomial::ChebyshevPolynomial<Complex> Poly(cheby_values, {-1, 1});
    std::cout<<exact/2-summand(0)-summ<<std::endl;
    std::cout<<Poly(0)<<std::endl;
    std::cout<<exact/2-summand(0)-Poly(1)<<" "<<Poly.ErrorEstimate()<<std::endl;
  }

}//namesapce
