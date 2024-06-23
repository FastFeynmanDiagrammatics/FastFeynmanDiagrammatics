#include<ffd.hpp>


using namespace ffd::sunset_propagator;
using namespace ffd::user_space;


std::size_t constexpr l_x = 6;
std::size_t constexpr l_y = 6;
std::size_t constexpr n_t = 60;
Real constexpr Beta = 5.;
Real constexpr mu0 =
  // -1.;
  // -0.7;
  -0.747691109658626;//beta5
// -.5;
  // -1;
  // -0.774574883247988;//beta20
  // -0.768519594956888;//beta10
Real constexpr t_hopping = 1.;
Real constexpr tp_hopping = -.3;
Real constexpr U = 5.6;
Real constexpr alpha_iter = 0.8;
Real constexpr alpha_iter_thir = 0.5;
std::size_t constexpr n_iter_bare =
  0;
std::size_t constexpr n_iter_bold =
  0;
std::size_t constexpr n_iter_ren =
  0;
std::size_t constexpr n_iter_third_bare =
  10;
std::size_t constexpr n_iter_third_bold =
  0;
std::size_t constexpr n_iter_third_ren =
  0;


int main(){
  using ffd::core_math::Pi;

  
  std::ostringstream SigmaR_br_name, SigmaR_bd_name,
    GR_bd_name, SigmaR_rs_name, GR_rs_name, SigmaR_th_name, GR_th_name, GR_br_name,
    common_part, end_part_br, end_part, name_br, name_bd, name_rs,  name_th, begin_part_G, begin_part_S;
  common_part<<std::setprecision(12)<<"_B"<<Beta;
  common_part<<"_U"<<U;
  common_part<<"_"<<(1<<l_x)<<"x"<<(1<<l_y);
  common_part<<"_Nc"<<n_t;
  common_part<<"_name";
  end_part<<std::setprecision(12)<<"_mu0"<<mu0;
  end_part<<"_tp"<<tp_hopping;
  end_part<<".dat";
  name_br<<"Bare2";
  name_bd<<"Bold2";
  name_rs<<"RSelf2";
  name_th<<"BareIII";
  begin_part_S<<"SigmaR";
  begin_part_G<<"GR____";
  SigmaR_br_name<<begin_part_S.str()<<common_part.str()<<name_br.str()<<end_part.str();
  SigmaR_bd_name<<begin_part_S.str()<<common_part.str()<<name_bd.str()<<end_part.str();
  SigmaR_rs_name<<begin_part_S.str()<<common_part.str()<<name_rs.str()<<end_part.str();
  SigmaR_th_name<<begin_part_S.str()<<common_part.str()<<name_th.str()<<end_part.str();
  GR_br_name<<begin_part_G.str()<<common_part.str()<<name_br.str()<<end_part.str();
  GR_bd_name<<begin_part_G.str()<<common_part.str()<<name_bd.str()<<end_part.str();
  GR_rs_name<<begin_part_G.str()<<common_part.str()<<name_rs.str()<<end_part.str();
  GR_th_name<<begin_part_G.str()<<common_part.str()<<name_th.str()<<end_part.str();
  std::ofstream SigmaR_br_file(SigmaR_br_name.str());
  std::ofstream SigmaR_bd_file(SigmaR_bd_name.str());
  std::ofstream SigmaR_rs_file(SigmaR_rs_name.str());
  std::ofstream SigmaR_th_file(SigmaR_th_name.str());
  std::ofstream GR_br_file(GR_br_name.str());
  std::ofstream GR_bd_file(GR_bd_name.str());
  std::ofstream GR_rs_file(GR_rs_name.str());
  std::ofstream GR_th_file(GR_th_name.str());

  
  auto const G0_kt =
    G0_k_time<l_x, l_y, n_t, ffd::l_array>(Beta, mu0, t_hopping, tp_hopping);
    
    
  auto const G0_rt = from_ktime_to_spacetime<l_x, l_y, n_t>(G0_kt);

  
  std::cerr<<"n(HF) = "<<-2.*G0_rt[0](Beta)<<"\n";

  
  auto [G_s_br, Sigma_s_br] =
    sunset_propagator_self_energy<l_x, l_y, n_t, n_iter_bare, ffd::phys::bare>(G0_rt,
									       U,
									       alpha_iter);
  auto G_s_br_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(G_s_br);
  auto Sigma_s_br_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(Sigma_s_br);
  std::cerr<<"n(semibare) = "<<-2.*G_s_br[0](Beta)<<std::endl;
  

  
  // auto Sigma_s_br_rt2 = sunset_self_energy<l_x, l_y, n_t>(G0_rt, U);
  // auto Sigma_s_br_kt2 = from_spacetime_to_ktime<l_x, l_y, n_t>(Sigma_s_br_rt2);
  // std::cerr<<"sigma "<<Sigma_s_br_kt2[0](0.)<<" "<<Sigma_s_br_kt[0](0.)<<"\n";
  
  auto [G_s_bd, Sigma_s_bd] =
    sunset_propagator_self_energy<l_x,
				  l_y,
				  n_t,
				  n_iter_bold,
				  ffd::phys::bold>(G0_rt,
						   U,
						   alpha_iter);
  auto G_s_bd_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(G_s_bd);
  auto Sigma_s_bd_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(Sigma_s_bd);
  std::cerr<<"n(semibold) = "<<-2.*G_s_bd[0](Beta)<<std::endl;


  auto [G_s_rs, Sigma_s_rs] = sunset_propagator_self_energy<l_x,
							    l_y,
							    n_t,
							    n_iter_ren,
							    ffd::phys::renormalized_self_energy>(G0_rt,
												 U,
												 alpha_iter);
  auto G_s_rs_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(G_s_rs);
  auto Sigma_s_rs_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(Sigma_s_rs);
  std::cerr<<"n(renormalized_self2) = "<<-2.*G_s_rs[0](Beta)<<std::endl;


  auto [G_s_th, Sigma_s_th] = third_order_propagator<l_x,
						     l_y,
						     n_t,
						     n_iter_third_bare,
						     ffd::phys::renormalized_self_energy>(G0_rt,
								      U,
								      alpha_iter);
  auto G_s_th_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(G_s_th);
  auto Sigma_s_th_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(Sigma_s_th);
  auto Sigma2_third_s = sunset_self_energy<l_x, l_y, n_t>(G_s_th, U);
  auto Sigma2_third_s_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(Sigma2_third_s);
  std::cerr<<"n(semibare_only3) = "<<-2.*G_s_th[0](Beta)<<std::endl;



  auto [G_s_TH, Sigma_s_TH] = third_order_propagator<l_x,
						     l_y,
						     n_t,
						     n_iter_third_bold,
						     ffd::phys::bold>(G0_rt,
								      U,
								      alpha_iter);
  auto G_s_TH_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(G_s_TH);
  auto Sigma_s_TH_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(Sigma_s_TH);
  auto Sigma2_THird_s = sunset_self_energy<l_x, l_y, n_t>(G_s_TH, U);
  auto Sigma2_THird_s_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(Sigma2_THird_s);
  std::cerr<<"n(semibold_only3) = "<<-2.*G_s_TH[0](Beta)<<std::endl;


  
  auto [G_s_TR, Sigma_s_TR] = third_order_propagator<l_x,
						     l_y,
						     n_t,
						     n_iter_third_ren,
						     ffd::phys::renormalized_self_energy>(G0_rt,
											  U,
											  alpha_iter);
  auto G_s_TR_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(G_s_TR);
  auto Sigma_s_TR_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(Sigma_s_TR);
  auto Sigma2_TR_s = sunset_self_energy<l_x, l_y, n_t>(G_s_TR, U);
  auto Sigma2_TR_s_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(Sigma2_TR_s);
  std::cerr<<"n(RSelf3) = "<<-2.*G_s_TR[0](Beta)<<std::endl;


  auto const chebyshev_nodes = ffd::chebyshev_polynomial_s::
    create_chebyshev_nodes<n_t>({0., Beta});
  for( std::size_t x=0; x < (1<<l_x); ++x){
    for( std::size_t y=0; y < (1<<l_y); ++y){
      std::size_t sp_ind = x*(1<<l_y)+y;
      for( std::size_t k=0; k < n_t; ++k){
	GR_br_file<<std::setprecision(12)<<G_s_br[sp_ind](chebyshev_nodes[k])<<'\n';
	GR_bd_file<<std::setprecision(12)<<G_s_bd[sp_ind](chebyshev_nodes[k])<<'\n';
	GR_rs_file<<std::setprecision(12)<<G_s_rs[sp_ind](chebyshev_nodes[k])<<'\n';
	GR_th_file<<std::setprecision(12)<<G_s_th[sp_ind](chebyshev_nodes[k])<<'\n';
	SigmaR_br_file<<std::setprecision(12)<<Sigma_s_br[sp_ind](chebyshev_nodes[k])<<'\n';
	SigmaR_bd_file<<std::setprecision(12)<<Sigma_s_bd[sp_ind](chebyshev_nodes[k])<<'\n';
	SigmaR_rs_file<<std::setprecision(12)<<Sigma_s_rs[sp_ind](chebyshev_nodes[k])<<'\n';
	SigmaR_th_file<<std::setprecision(12)<<Sigma_s_th[sp_ind](chebyshev_nodes[k])<<'\n';
      }
    }
  }
  
  
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

  
  std::ofstream spectral_omega0_HF("spectral_omega0_HF.dat");
  std::ofstream spectral_omega01_HF("spectral_omega01_HF.dat");
  std::ofstream spectral_omega0_bare("spectral_omega0_bare.dat");
  std::ofstream spectral_omega01_bare("spectral_omega01_bare.dat");
  std::ofstream spectral_omega0_RSelf("spectral_omega0_RSelf.dat");
  std::ofstream spectral_omega0_BareThird("spectral_omega0_BareThird.dat");
  std::ofstream spectral_omega01_RSelf("spectral_omega01_RSelf.dat");
  std::ofstream spectral_omega0_bold("spectral_omega0_bold.dat");
  std::ofstream spectral_omega01_bold("spectral_omega01_bold.dat");
  std::ofstream imagpart_omega0_HF("imagpart_omega0_HF.dat");
  std::ofstream imagpart_omega01_HF("imagpart_omega01_HF.dat");
  std::ofstream imagpart_omega0_bare("imagpart_omega0_bare.dat");
  std::ofstream imagpart_omega01_bare("imagpart_omega01_bare.dat");
  std::ofstream imagpart_omega0_bold("imagpart_omega0_bold.dat");
  std::ofstream imagpart_omega01_bold("imagpart_omega01_bold.dat");
  std::ofstream imagpart_omega0_RSelf("imagpart_omega0_RSelf.dat");
  std::ofstream imagpart_omega01_RSelf("imagpart_omega01_RSelf.dat");
  std::ofstream imagpart_omega0_BareThird("imagpart_omega0_BareThird.dat");
  std::ofstream imagpart_omega01_BareThird("imagpart_omega01_BareThird.dat");
  std::ofstream imagpart_omega0_BareThirdSecondOrder("imagpart_omega0_BareThirdSecondOrder.dat");


  
  std::ofstream realpart_omega0_HF("realpart_omega0_HF.dat");
  std::ofstream realpart_omega01_HF("realpart_omega01_HF.dat");
  std::ofstream realpart_omega0_bare("realpart_omega0_bare.dat");
  std::ofstream realpart_omega01_bare("realpart_omega01_bare.dat");
  std::ofstream realpart_omega0_bold("realpart_omega0_bold.dat");
  std::ofstream realpart_omega01_bold("realpart_omega01_bold.dat");
  std::ofstream realpart_omega0_RSelf("realpart_omega0_RSelf.dat");
  std::ofstream realpart_omega01_RSelf("realpart_omega01_RSelf.dat");
  std::ofstream realpart_omega0_BareThird("realpart_omega0_BareThird.dat");
  std::ofstream realpart_omega01_BareThird("realpart_omega01_BareThird.dat");


  std::ofstream fermisurface_omega0_HF("fermisurface_omega0_HF.dat");
  std::ofstream fermisurface_omega01_HF("fermisurface_omega01_HF.dat");
  std::ofstream fermisurface_omega0_bare("fermisurface_omega0_bare.dat");
  std::ofstream fermisurface_omega01_bare("fermisurface_omega01_bare.dat");
  std::ofstream fermisurface_omega0_bold("fermisurface_omega0_bold.dat");
  std::ofstream fermisurface_omega01_bold("fermisurface_omega01_bold.dat");
  std::ofstream fermisurface_omega0_RSelf("fermisurface_omega0_RSelf.dat");
  std::ofstream fermisurface_omega01_RSelf("fermisurface_omega01_RSelf.dat");
  std::ofstream fermisurface_omega0_BareThird("fermisurface_omega0_BareThird.dat");
  std::ofstream fermisurface_omega01_BareThird("fermisurface_omega01_BareThird.dat");

  
  
  Real err_max1 = 0., err_max2 = 0., err_max3 = 0.;
  std::array<ffd::l_array<Complex, size(G0_kt)>, 2> Sigma_s_br_k;
  std::array<ffd::l_array<Complex, size(G0_kt)>, 2> Sigma_s_bd_k;
  std::array<ffd::l_array<Complex, size(G0_kt)>, 2> Sigma_s_rs_k;
  std::array<ffd::l_array<Complex, size(G0_kt)>, 2> Sigma_s_th_k;
  std::array<ffd::l_array<Complex, size(G0_kt)>, 2> Sigma2_third_s_k;
  std::array<ffd::l_array<Complex, size(G0_kt)>, 2> G_s_bd_k;
  std::array<ffd::l_array<Complex, size(G0_kt)>, 2> G_s_rs_k;
  std::array<ffd::l_array<Complex, size(G0_kt)>, 2> G_s_th_k;
  std::array<ffd::l_array<Complex, size(G0_kt)>, 2> G_s_br_k;
  std::array<ffd::l_array<Complex, size(G0_kt)>, 2> G0_k;
  for(std::size_t omega: {0, 1}){
    for(std::size_t k=0; k<size(G0_kt); ++k){
      auto [kx, ky] = ffd::math_tools::ModDiv(k, 1ul<<l_x);
      using ffd::core_math::Pi;
      Real Kx = 2.*kx*ffd::core_math::Pi/(1ul<<l_x);
      // if(Kx>Pi){
      // 	Kx -= 2.*Pi;
      // }
      Real Ky = 2.*ky*ffd::core_math::Pi/(1ul<<l_y);
      // if(Ky>Pi){
      // 	Ky -= 2.*Pi;
      // }

      
      Sigma_s_br_k[omega][k] = 0.;
      Sigma_s_bd_k[omega][k] = 0.;
      Sigma_s_rs_k[omega][k] = 0.;
      Sigma_s_th_k[omega][k] = 0.;
      Sigma2_third_s_k[omega][k] = 0.;
      G_s_br_k[omega][k] = 0.;
      G_s_bd_k[omega][k] = 0.;
      G_s_rs_k[omega][k] = 0.;
      G_s_th_k[omega][k] = 0.;
      G0_k[omega][k] = 0.;
      for(std::size_t j=0; j<n_t; ++j){
  	Sigma_s_br_k[omega][k] += Beta*Sigma_s_br_kt[k].coef[j]*matsu_coef[omega][j];
	Sigma_s_bd_k[omega][k] += Beta*Sigma_s_bd_kt[k].coef[j]*matsu_coef[omega][j];
	Sigma_s_rs_k[omega][k] += Beta*Sigma_s_rs_kt[k].coef[j]*matsu_coef[omega][j];
	Sigma_s_th_k[omega][k] += Beta*Sigma_s_th_kt[k].coef[j]*matsu_coef[omega][j];
	Sigma2_third_s_k[omega][k] += Beta*Sigma2_third_s_kt[k].coef[j]*matsu_coef[omega][j];
	G_s_bd_k[omega][k] += Beta*G_s_bd_kt[k].coef[j]*matsu_coef[omega][j];
	G_s_rs_k[omega][k] += Beta*G_s_rs_kt[k].coef[j]*matsu_coef[omega][j];
	G_s_th_k[omega][k] += Beta*G_s_th_kt[k].coef[j]*matsu_coef[omega][j];
	G_s_br_k[omega][k] += Beta*G_s_br_kt[k].coef[j]*matsu_coef[omega][j];
	G0_k[omega][k] += Beta*G0_kt[k].coef[j]*matsu_coef[omega][j];
      }
      
      // if(omega==1){
      // 	std::cerr<<k<<" 0= "<<Sigma_s_rs_k[omega][k]<<std::endl;
      // }
      
      using std::pow, std::abs;
      if(omega==0 && kx <= (1ul<<(l_x-1)) && ky <= (1ul<<(l_y-1)) ){
	spectral_omega0_HF<<Kx<<" "<<Ky<<" "<<imag(G0_k[0][k])<<"\n";
	spectral_omega0_bare<<Kx<<" "<<Ky<<" "<<
	  -imag(Sigma_s_br_k[0][k])/
	  (pow(real(1./G0_k[0][k]-Sigma_s_br_k[0][k]), 2)+pow(imag(Sigma_s_br_k[0][k]),2))<<"\n";
	spectral_omega0_bold<<Kx<<" "<<Ky<<" "<<
	  -imag(Sigma_s_bd_k[0][k])/
	  (pow(real(1./G0_k[0][k]-Sigma_s_bd_k[0][k]), 2)+pow(imag(Sigma_s_bd_k[0][k]),2))<<"\n";
	spectral_omega0_RSelf<<Kx<<" "<<Ky<<" "<<
	  -imag(Sigma_s_rs_k[0][k])/
	  (pow(real(1./G0_k[0][k]-Sigma_s_rs_k[0][k]), 2)+pow(imag(Sigma_s_rs_k[0][k]),2))<<"\n";
	spectral_omega0_BareThird<<Kx<<" "<<Ky<<" "<<
	  -imag(Sigma_s_th_k[0][k])/
	  (pow(real(1./G0_k[0][k]-Sigma_s_th_k[0][k]), 2)+pow(imag(Sigma_s_th_k[0][k]),2))<<"\n";

	
	imagpart_omega0_bare<<Kx<<" "<<Ky<<" "<<
	    imag(Sigma_s_br_k[0][k])<<"\n";
	imagpart_omega0_bold<<Kx<<" "<<Ky<<" "<<
	    imag(Sigma_s_bd_k[0][k])<<"\n";
	imagpart_omega0_RSelf<<Kx<<" "<<Ky<<" "<<
	    imag(Sigma_s_rs_k[0][k])<<"\n";
	imagpart_omega0_BareThird<<Kx<<" "<<Ky<<" "<<
	    imag(Sigma_s_th_k[0][k])<<"\n";
	imagpart_omega0_BareThirdSecondOrder<<Kx<<" "<<Ky<<" "<<
	    imag(Sigma2_third_s_k[0][k])<<"\n";
	
	  
	fermisurface_omega0_bare<<Kx<<" "<<Ky<<" "<<
	  -real(1./G0_k[0][k]-Sigma_s_br_k[0][k])<<"\n";
	fermisurface_omega0_HF<<Kx<<" "<<Ky<<" "<<
	  -real(1./G0_k[0][k])<<"\n";
	fermisurface_omega0_bold<<Kx<<" "<<Ky<<" "<<
	  -real(1./G0_k[0][k]-Sigma_s_bd_k[0][k])<<"\n";
	fermisurface_omega0_RSelf<<Kx<<" "<<Ky<<" "<<
	  -real(1./G0_k[0][k]-Sigma_s_rs_k[0][k])<<"\n";

	

	realpart_omega0_bare<<Kx<<" "<<Ky<<" "<<
	  -real(-Sigma_s_br_k[0][k])<<"\n";
	realpart_omega0_HF<<Kx<<" "<<Ky<<" "<<
	  -std::real(0.)<<"\n";
	realpart_omega0_bold<<Kx<<" "<<Ky<<" "<<
	  -real(-Sigma_s_bd_k[0][k])<<"\n";
	realpart_omega0_RSelf<<Kx<<" "<<Ky<<" "<<
	  -real(-Sigma_s_rs_k[0][k])<<"\n";

	}
	
      auto g0 = 1./G0_k[omega][k];
      auto g1 = std::abs(1./G_s_br_k[omega][k]+Sigma_s_br_k[omega][k]-g0);
      auto g2 = std::abs(1./G_s_bd_k[omega][k]+Sigma_s_bd_k[omega][k]-g0);
      auto g3 = std::abs(1./G_s_rs_k[omega][k]+Sigma_s_rs_k[omega][k]-g0);
      // std::cerr<<g0<<" "<<
      // 	g1<<" "<<g2<<'	
      if(g1>err_max1){
	err_max1 = g1;
      }
      if(g2>err_max2){
	err_max2= g2;
      }
      if(g3>err_max3){
	err_max3= g3;
      }
	    
    }
  }

  std::cerr<<"err_max_bare = "<<err_max1<<" err_max_bold="<<err_max2<<" err_max_RSelf="<<err_max3<<"\n";


  // for(std::size_t kx=0; kx<=(1<<(l_x-1)); ++kx){
  //   Real const Kx = 2.*kx*ffd::core_math::Pi/(1ul<<l_x);
  //   auto sigma_0 = Sigma_omega_k[0][kx+kx*(1<<l_x)];
  //   auto sigma_1 = Sigma_omega_k[1][kx+kx*(1<<l_x)];
  //   std::cerr<<std::scientific<<Kx<<"\t"<<
  //     real(sigma_0)<<"\t"<<
  //     imag(sigma_0)<<"\t"<<
  //     imag(sigma_1)<<"\n";
    
  // }

  
  // std::ofstream SigmaI_omega0("SigmaI_omega0.dat");
  // for(std::size_t kx=0; kx<=(1<<(l_x-1)); ++kx){
  //   Real const Kx = 2.*kx*ffd::core_math::Pi/(1ul<<l_x);
    //   for(std::size_t ky=0; ky<=(1<<(l_y-1)); ++ky){
  //     Real const Ky = 2.*ky*ffd::core_math::Pi/(1ul<<l_y);
  //     auto sigma_0 = Sigma_omega_k[0][kx+ky*(1<<l_x)];
  //     SigmaI_omega0<<Kx<<" "<<Ky<<" "<<imag(sigma_0)<<"\n";
  //   }
  // }
  

  
  // std::ofstream spectral_omega0("spectral_omega0.dat");
  // for(std::size_t kx=0; kx<=(1<<(l_x-1)); ++kx){
  //   Real const Kx = 2.*kx*ffd::core_math::Pi/(1ul<<l_x);
  //   for(std::size_t ky=0; ky<=(1<<(l_y-1)); ++ky){
  //     Real const Ky = 2.*ky*ffd::core_math::Pi/(1ul<<l_y);
  //     Complex const Iomega = Complex(0., (2*0+1)*Pi/Beta);
  //     Real const cos_kx = cos(Kx);
  //     Real const cos_ky = cos(Ky);
  //     Real const xi_k = -2*t_hopping*(cos_kx+cos_ky) -4*tp_hopping*cos_kx*cos_ky - mu0;
  //     Complex const sigma_k_omega0 = Sigma_omega_k[0][kx+ky*(1<<l_x)];
  //     Complex const G_k_omega0 = 1./(Iomega-xi_k-sigma_k_omega0);
  //     spectral_omega0<<Kx<<" "<<Ky<<" "<<imag(G_k_omega0)<<"\n";
  //   }
  // }
  
  
  
}
