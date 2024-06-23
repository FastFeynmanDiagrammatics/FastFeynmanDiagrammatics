#include <ffd.hpp>

using namespace ffd::user_space;

Real constexpr Beta = 1.;
Real constexpr Mu = -0.4054651081;
  // -0.31128232741126193872;
Real constexpr U = 1.;
int constexpr omega_min = 0;
int constexpr omega_max = 0;
int constexpr order = 13;
std::size_t constexpr num_iter = 1ul<<27;
bool constexpr no_tadpoles = true;


int main(){
  using ffd::core_math::Pi;
  auto H = [=](auto n, auto mu, auto z){
	     return -mu*Real(n[0]+n[1]) + z*Real(n[0]*n[1]);
	   };

  
  auto Boltz = [=](auto n, auto mu, auto z){
		 using std::exp;
		 return exp(-Beta*H(n, mu, z));
	       };

  
  auto Z = [=](auto mu, auto z){
	     using std::exp;
	     decltype(z) ret = 0.;
	     for(int n_up: {0, 1}){
	       for(int n_do: {0, 1}){
		 ret += Boltz(std::array<int, 2>{n_up, n_do}, mu, z);
	       }
	     }
	     return ret;
	   };

  
  auto z1 = [=](auto mu, auto z){
	     decltype(z) ret = 0.;
	     for(int n_up: {0, 1}){
	       for(int n_do: {0, 1}){
		 ret += Real(n_up)*Boltz(std::array<int, 2>{n_up, n_do}, mu, z);
	       }
	     }
	     return ret/Z(mu, z);
	   };
 
  
  Real const magic_alpha = no_tadpoles ? z1(Mu, 0.) : 0.;
  std::cerr<< "n0 = " << 2.*magic_alpha << std::endl;

  
  
  auto G0 = [=](auto tau){
	      Real sign = -1;
	      while(tau < 0){
		tau += Beta;
		sign = -sign;
	      }
	      while(tau > Beta){
		tau -= Beta;
		sign = -sign;
	      }
	      return sign*(1-magic_alpha)*std::exp(tau*Mu);
	    };
  
  auto G0_omega = [=](int omega){
		    return 1./(Complex(0., (.5+omega)*2*Pi/Beta)+Mu);
		  };

  auto Sigma = [=](auto mu, auto z, int n_omega){
		 auto I_omega = Complex(0., (n_omega+.5)*2.*Pi/Beta);
		 auto Z1 = z1(mu, z);
		 return Z1*(1.-Z1)*z*z/(I_omega+mu-(1.-Z1)*z);
	       };
  
  
  // for(int j=omega_min; j<=omega_max; ++j){
  //   Complex g0 = 0.;
  //   for(std::size_t t_MC=0; t_MC<num_iter; ++t_MC){
  //     Real tau = Proba()*Beta;
  //     g0 += std::exp(Complex(0., (0.5+j)*2.*Pi*tau/Beta))*G0(tau);
  //   }
  //   g0 /= Real(num_iter);
  //   std::cerr<<"G^{(0)}_{IOmega_"<<j<<"} = "<<g0<<std::endl;
  //   std::cerr<<"exact G^{(0)}_{IOmega_"<<j<<"} = "<<G0_omega(j)<<std::endl;
  // }


  // for(int j=omega_min; j<=omega_max; ++j){
  //   Complex Xi_2 = 0.;
  //   Real omega = (.5+j)*2.*Pi/Beta;
  //   for(std::size_t t_MC=0; t_MC<num_iter; ++t_MC){
  //     Real tau = Proba()*Beta;
  //     Xi_2 += Beta*std::exp(Complex(0., (0.5+j)*2.*Pi*tau/Beta))*G0(tau)*G0(tau)*G0(-tau);
  //   }
  //   Xi_2 /= Real(num_iter);
  //   std::cerr<<"Xi^{(2)}_{IOmega_"<<j<<"} = "<<Xi_2<<std::endl;
  //   std::cerr<<"exact Xi^{(2)}_{IOmega_"<<j<<"} = "<<magic_alpha/(Complex(0., omega)+Mu)<<std::endl;
  // }


  for(int j=omega_min; j<=omega_max; ++j){
    auto sigma_der = ffd::contour_derivative::
      ContourDerivative([=](Complex z){return Sigma(Mu+magic_alpha*z, z, j);},
			1.);
    std::cerr<<"Sigma_{Iomega_"<<j<<"} = "<<Sigma(Mu+magic_alpha*U, U, j)<<std::endl;
    for(int k=2; k<=order; ++k){
      Real fac = ffd::core_math::Factorial(k);
            std::cerr<<"Re Sigma^{("<<k<<")}_{Iomega_"<<j<<"} "<<std::real(sigma_der(k)/fac)<<"\n";
    }
    std::cerr<<"\n";
      for(int k=2; k<=order; ++k){
      Real fac = ffd::core_math::Factorial(k);

      std::cerr<<"Im Sigma^{("<<k<<")}_{Iomega_"<<j<<"} "<<std::imag(sigma_der(k)/fac)<<"\n";
    }

    std::cerr<<"\n";
  }

  
    
  
  
}
