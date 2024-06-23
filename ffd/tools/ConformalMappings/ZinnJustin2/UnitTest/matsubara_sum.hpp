

namespace ffd::conformal_mappings::unit_test{

  void matsubara_sum(){
    using ffd::core_math::Pi;
    Complex I = Complex(0., 1.);
    const int max_order = 50;
    const Real Beta = 1;
    const Real xi = 1;
    //    const Real tau = 1;
    std::vector<Complex> coef_z(max_order);
    auto omega_n = [Beta](int n) -> Real {return 2*Pi*(n+.5)/Beta;};
    //    auto coef_gen = [](int n) -> Complex {return 1.;};//(I*omega_n(n)*omega_n(n) - xi);};
    auto coef_gen = [I, omega_n, xi,Beta](int n) -> Complex
		    {return (1./Beta)*(1./(I*omega_n(n) - xi)+I/omega_n(n));};
    Complex res = 0.;
    Real tau = Beta/2;
    Complex z0 = exp(-2.*I*Pi*tau/Beta);
    Complex sqrt_z0 = std::sqrt(z0);
    for(int j=0; j < max_order; ++j){
      coef_z[j] = coef_gen(j)*std::pow(-1., j);
      res += coef_gen(j)*std::pow(z0, j);
    }
    ZinnJustin2<Complex> ZJ(coef_z);
    std::cerr<<2*std::real(sqrt_z0*ZJ(-z0))<<" "<<2*std::real(sqrt_z0*res)<<std::endl;
    std::cerr<<(-exp(-tau*xi)/(1+exp(-Beta*xi))+.5)<<std::endl;

    std::vector<Real> coef_z2(max_order, 0);
    for(int j=1; j < max_order; ++j){
      coef_z2[j] = pow(-1, j+1)/(j);
    }
    ZinnJustin2<Real> ZJ2(coef_z2);
    Complex z0_log = exp(I*2.);
    std::cout<<ZJ2(z0_log)/log(1.+z0_log)<<std::endl;
  }

}//namespace
