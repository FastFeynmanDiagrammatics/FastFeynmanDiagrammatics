#include<ffd.hpp>


using namespace ffd::sunset_propagator;
using namespace ffd::user_space;

std::size_t constexpr l_x = 8;
std::size_t constexpr l_y = 8;
std::size_t constexpr n_t = 50;
Real constexpr Beta = 5.;
Real constexpr mu0 =-0.747691109658626;//beta5
  //-0.774574883247988;//beta20
  //-0.768519594956888;//beta10
Real constexpr t_hopping = 1.;
Real constexpr tp_hopping = -.3;
Real constexpr U = 5.6;


int main(){
  using ffd::core_math::Pi;
  
  auto const G0_kt =
    G0_k_time<l_x, l_y, n_t, ffd::l_array>(Beta, mu0, t_hopping, tp_hopping);
    
    
  auto const G0_rt = from_ktime_to_spacetime<l_x, l_y, n_t>(G0_kt);

  std::cerr<<"n0 = "<<-2.*G0_rt[0](Beta)<<"\n";


  auto const Sigma_rt = sunset_self_energy<l_x, l_y, n_t>(G0_rt, U);


  auto const Sigma_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(Sigma_rt);


  std::array<std::array<Complex, n_t>, 2> matsu_coef;
  for(std::size_t omega: {0, 1}){
    matsu_coef[omega] = ffd::chebyshev_polynomial_s::matsubara_transform<n_t>(omega);
  }


  // ffd::l_array<Complex, size(Sigma_kt)> G0_ko;
  // for(std::size_t k=0; k<size(Sigma_kt); ++k){
  //   G0_ko[k] = 0.;
  //   for(std::size_t j=0; j<n_t; ++j){
  //     G0_ko[k] += Beta*G0_kt[k].coef[j]*matsu_coef[j];
  //   }
  //   Complex const Iomega = Complex(0., (2*0+1)*Pi/Beta);
  //   auto [kx, ky] = ffd::math_tools::ModDiv(k, (1<<l_x));
  //   Real const cos_kx = cos(kx*2.*Pi/(1<<l_x));
  //   Real const cos_ky = cos(ky*2.*Pi/(1<<l_y));
  //   Real const xi_k = -2*t_hopping*(cos_kx+cos_ky) -4*tp_hopping*cos_kx*cos_ky - mu0;
  //   std::cerr<<G0_ko[k]<<" "<<1./(Iomega-xi_k)<<std::endl;
  // }
  

  
  std::array<ffd::l_array<Complex, size(Sigma_kt)>, 2> Sigma_omega_k;
  for(std::size_t omega: {0, 1}){
    for(std::size_t k=0; k<size(Sigma_kt); ++k){
      Sigma_omega_k[omega][k] = 0.;
      for(std::size_t j=0; j<n_t; ++j){
	Sigma_omega_k[omega][k] += Beta*Sigma_kt[k].coef[j]*matsu_coef[omega][j];
      }
    }
  }
  


  for(std::size_t kx=0; kx<=(1<<(l_x-1)); ++kx){
    Real const Kx = 2.*kx*ffd::core_math::Pi/(1ul<<l_x);
    auto sigma_0 = Sigma_omega_k[0][kx+kx*(1<<l_x)];
    auto sigma_1 = Sigma_omega_k[1][kx+kx*(1<<l_x)];
    std::cerr<<std::scientific<<Kx<<"\t"<<
      real(sigma_0)<<"\t"<<
      imag(sigma_0)<<"\t"<<
      imag(sigma_1)<<"\n";
    
  }

  
  std::ofstream SigmaI_omega0("SigmaI_omega0.dat");
  for(std::size_t kx=0; kx<=(1<<(l_x-1)); ++kx){
    Real const Kx = 2.*kx*ffd::core_math::Pi/(1ul<<l_x);
    for(std::size_t ky=0; ky<=(1<<(l_y-1)); ++ky){
      Real const Ky = 2.*ky*ffd::core_math::Pi/(1ul<<l_y);
      auto sigma_0 = Sigma_omega_k[0][kx+ky*(1<<l_x)];
      SigmaI_omega0<<Kx<<" "<<Ky<<" "<<imag(sigma_0)<<"\n";
    }
  }
  

  
  std::ofstream spectral_omega0("spectral_omega0.dat");
  for(std::size_t kx=0; kx<=(1<<(l_x-1)); ++kx){
    Real const Kx = 2.*kx*ffd::core_math::Pi/(1ul<<l_x);
    for(std::size_t ky=0; ky<=(1<<(l_y-1)); ++ky){
      Real const Ky = 2.*ky*ffd::core_math::Pi/(1ul<<l_y);
      Complex const Iomega = Complex(0., (2*0+1)*Pi/Beta);
      Real const cos_kx = cos(Kx);
      Real const cos_ky = cos(Ky);
      Real const xi_k = -2*t_hopping*(cos_kx+cos_ky) -4*tp_hopping*cos_kx*cos_ky - mu0;
      Complex const sigma_k_omega0 = Sigma_omega_k[0][kx+ky*(1<<l_x)];
      Complex const G_k_omega0 = 1./(Iomega-xi_k-sigma_k_omega0);
      spectral_omega0<<Kx<<" "<<Ky<<" "<<imag(G_k_omega0)<<"\n";
    }
  }
  
  
  
}
