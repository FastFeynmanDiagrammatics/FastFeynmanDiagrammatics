#include<ffd.hpp>
#include"parameters.hpp"

using namespace ffd::user_space;

using namespace ffd::sunset_propagator;
using namespace ffd::user_space::parameters;


std::size_t constexpr size_x_y= (1ul<<(l_x+l_y));
std::size_t constexpr size_x = 1ul<<l_x;
std::size_t constexpr size_y = 1ul<<l_y;
 
using std::cos, std::exp, std::sin, std::cosh;
using ffd::core_math::Pi;


#include"kt_propagator.hpp"


int main(){
 
  #include"lambdas.hpp"
  
  // int ch_flat_print = 100;
  // for(int j=0; j<ch_flat_print; ++j){
  //   Real tau = (j+.5)*Beta/ch_flat_print;
  //   std::cerr<<tau<<' '<<flat2_ch_t(tau)<<'\n';
  // }
  

  
  
  fermi_G_t<l_x, l_y, n_t> G0_kt, G0_rt;
  Real n0 = 0.25;
  Real mu_0 = mu - U*(n0+alpha_shift);
  for(int j=0; j<200; ++j){
    mu_0 = .5*mu_0 + .5*(mu - U*(n0+alpha_shift));
    for(std::size_t nx=0; nx<size_x; ++nx){
      for(std::size_t ny=0; ny<size_y; ++ny){
	Real const kx = 2.*nx*Pi/size_x;
	Real const ky = 2.*ny*Pi/size_y;
	std::size_t const index = ny+size_y*nx;
	G0_kt[index] = fermi_cheby_s<n_t>([=](Real t){return green_function0_kt_mu0(kx, ky, t, mu_0);},
					  {0., Beta});
      }
    }
    G0_rt = from_ktime_to_spacetime<l_x, l_y, n_t>(G0_kt);
    n0 =   -G0_rt[0](Beta);
    // std::cerr<<2.*n0<<'\n';
  }
  
  
  fermi_G_t<l_x, l_y, n_t> YRZ_si_kt, BCS_si_kt, BCD_si_kt, mix_si_kt, FLAT_si_kt;
  for(std::size_t nx=0; nx<size_x; ++nx){
    for(std::size_t ny=0; ny<size_y; ++ny){
      Real const kx = 2.*nx*Pi/size_x;
      Real const ky = 2.*ny*Pi/size_y;
      std::size_t const index = ny+size_y*nx;
      YRZ_si_kt[index] = fermi_cheby_s<n_t>([=](Real t){return YRZsi_kt(kx, ky, t);},
					    {0., Beta});
      BCS_si_kt[index] = fermi_cheby_s<n_t>([=](Real t){return BCSsi_kt_mu0(kx, ky, t, mu_0);},
					    {0., Beta});
      BCD_si_kt[index] = fermi_cheby_s<n_t>([=](Real t){return BCDsi_kt_mu0(kx, ky, t, mu_0);},
					    {0., Beta});
      mix_si_kt[index] = fermi_cheby_s<n_t>([=](Real t){return
	    BCSsi_kt_mu0(kx, ky, t, mu_0)+FLATsi_kt_mu0(kx, ky, t, mu_0);}
	    // mix_BCS*BCSsi_kt_mu0(kx, ky, t, mu_0)+
	    // mix_BCD*BCDsi_kt_mu0(kx, ky, t, mu_0)+mix_YRZ*YRZsi_kt(kx, ky, t);},
	,{0., Beta});
      FLAT_si_kt[index] = fermi_cheby_s<n_t>([=](Real t){return FLAT2si_kt_mu0(kx, ky, t, mu_0)
	    + BCSsi_kt_mu0(kx, ky, t, mu_0);},
					    {0., Beta});
    }
  }
  std::array<Real, size_x> G0_norm_nonlocal;
  G0_norm_nonlocal.fill(0.);
  for(std::size_t x=0; x<size_x; ++x){
    for(std::size_t y=0; y<size_y; ++y){
      std::size_t const norm_max_r = std::max(x, y);
      G0_norm_nonlocal[norm_max_r] += norm_max(G0_rt[y+size_y*x]);
    }
  }
  std::cerr<<"n_0 = "<<-2.*G0_rt[0](Beta)<<" norm_non_local= ";
  int print_norm_values = 10;
  for(int j=0; j<print_norm_values; ++j){
    std::cerr<<std::setprecision(2)<<G0_norm_nonlocal[j]<<' ';
  }
  std::cerr<<'\n'<<std::setprecision(5);


  
  auto const YRZ_si_rt = from_ktime_to_spacetime<l_x, l_y, n_t>(YRZ_si_kt);
  auto YRZ_G_rt = G0_rt;
  auto const BCS_si_rt = from_ktime_to_spacetime<l_x, l_y, n_t>(BCS_si_kt);
  auto BCS_G_rt = G0_rt;
  auto const BCD_si_rt = from_ktime_to_spacetime<l_x, l_y, n_t>(BCD_si_kt);
  auto BCD_G_rt = G0_rt;
  auto mix_si_rt = from_ktime_to_spacetime<l_x, l_y, n_t>(mix_si_kt);
  auto mix_G_rt = G0_rt;
  auto const FLAT_si_rt = from_ktime_to_spacetime<l_x, l_y, n_t>(FLAT_si_kt);
  auto FLAT_G_rt = G0_rt;
  
  
  for(std::size_t iter=0; iter<iter_YRZ; ++iter){
    YRZ_G_rt = iteration_dyson_equation<l_x, l_y, n_t>(G0_rt, YRZ_si_rt, YRZ_G_rt, a_iter_YRZ);
  }
  Real YRZ_norm_nonlocal = 0.;
  for(std::size_t r=1; r<size_x_y; ++r){
    YRZ_norm_nonlocal += -2.*std::abs(YRZ_G_rt[r](Beta));
  }
  if(iter_YRZ>0){
    std::cerr<<"n_YRZ = "<<-2.*YRZ_G_rt[0](Beta)
	     <<" norm_nonlocal = "<<YRZ_norm_nonlocal
	     <<" iter = "<<iter_YRZ<<'\n';
  }
  auto const YRZ_si2_rt = sunset_self_energy<l_x, l_y, n_t>(YRZ_G_rt, U);
  auto const YRZ_si2_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(YRZ_si2_rt);
  auto const YRZ_si3_rt = third_order_self_energy<l_x, l_y, n_t>(YRZ_G_rt, U);
  auto const YRZ_si3_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(YRZ_si3_rt);


  for(std::size_t iter=0; iter<iter_BCS; ++iter){
    BCS_G_rt = iteration_dyson_equation<l_x, l_y, n_t>(G0_rt, BCS_si_rt, BCS_G_rt, a_iter_BCS);
  }
  std::array<Real, size_x> BCS_norm_nonlocal;
  BCS_norm_nonlocal.fill(0.);
  for(std::size_t x=0; x<size_x; ++x){
    for(std::size_t y=0; y<size_y; ++y){
      std::size_t const norm_max_r = std::max(x, y);
      BCS_norm_nonlocal[norm_max_r] += norm_max(BCS_G_rt[y+size_y*x]);
    }
  }
  std::cerr<<"n_BCS = "<<-2.*BCS_G_rt[0](Beta)<<" norm_non_local= ";
  for(int j=0; j<print_norm_values; ++j){
    std::cerr<<std::setprecision(2)<<BCS_norm_nonlocal[j]<<' ';
  }
  std::cerr<<'\n'<<std::setprecision(5);


  auto const BCS_si2_rt = sunset_self_energy<l_x, l_y, n_t>(BCS_G_rt, U);
  auto const BCS_si2_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(BCS_si2_rt);
  auto const BCS_si3_rt = third_order_self_energy<l_x, l_y, n_t>(BCS_G_rt, U);
  auto const BCS_si3_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(BCS_si3_rt);


  for(std::size_t iter=0; iter<iter_BCD; ++iter){
    BCD_G_rt = iteration_dyson_equation<l_x, l_y, n_t>(G0_rt, BCD_si_rt, BCD_G_rt, a_iter_BCD);
  }
  if(iter_BCD>0){
    std::cerr<<"n_BCD = "<<-2.*BCD_G_rt[0](Beta)<<" iter = "<<iter_BCD<<'\n';
  }
  auto const BCD_si2_rt = sunset_self_energy<l_x, l_y, n_t>(BCD_G_rt, U);
  auto const BCD_si2_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(BCD_si2_rt);


  
  for(std::size_t iter=0; iter<iter_FLAT; ++iter){
    FLAT_G_rt = iteration_dyson_equation<l_x, l_y, n_t>(G0_rt, FLAT_si_rt, FLAT_G_rt, a_iter_FLAT);
  }
  auto const FLAT_si2_rt = sunset_self_energy<l_x, l_y, n_t>(FLAT_G_rt, U);
  auto const FLAT_si2_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(FLAT_si2_rt);
  auto const FLAT_si3_rt = third_order_self_energy<l_x, l_y, n_t>(FLAT_G_rt, U);
  auto const FLAT_si3_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(FLAT_si3_rt);


  std::array<Real, size_x> FLAT_norm_nonlocal;
  FLAT_norm_nonlocal.fill(0.);
  for(std::size_t x=0; x<size_x; ++x){
    for(std::size_t y=0; y<size_y; ++y){
      std::size_t const norm_max_r = std::max(x, y);
      FLAT_norm_nonlocal[norm_max_r] += norm_max(FLAT_G_rt[y+size_y*x]);
    }
  }
  std::cerr<<"n_FLAT = "<<-2.*FLAT_G_rt[0](Beta)<<" norm_non_local= ";
  for(int j=0; j<print_norm_values; ++j){
    std::cerr<<std::setprecision(2)<<FLAT_norm_nonlocal[j]<<' ';
  }
  std::cerr<<'\n'<<std::setprecision(5);



  
  mix_G_rt = FLAT_G_rt;
  for(std::size_t iter=0; iter<iter_mix; ++iter){
    mix_si_rt = third_order_self_energy<l_x, l_y, n_t>(mix_G_rt, U);
    mix_G_rt = iteration_dyson_equation<l_x, l_y, n_t>(FLAT_G_rt, mix_si_rt, mix_G_rt, a_iter_mix);
  }
  if(iter_mix>0){
      std::array<Real, size_x> mix_norm_nonlocal;
  mix_norm_nonlocal.fill(0.);
  for(std::size_t x=0; x<size_x; ++x){
    for(std::size_t y=0; y<size_y; ++y){
      std::size_t const norm_max_r = std::max(x, y);
      mix_norm_nonlocal[norm_max_r] += norm_max(mix_G_rt[y+size_y*x]);
    }
  }
  std::cerr<<"n_mix = "<<-2.*mix_G_rt[0](Beta)<<" norm_non_local= ";
  for(int j=0; j<print_norm_values; ++j){
    std::cerr<<std::setprecision(2)<<mix_norm_nonlocal[j]<<' ';
  }
  std::cerr<<'\n'<<std::setprecision(5);

  }
  auto const mix_si2_rt = sunset_self_energy<l_x, l_y, n_t>(mix_G_rt, U);
  auto const mix_si2_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(mix_si2_rt);
  auto const mix_si3_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(mix_si_rt);


  
  auto ch_to_mts = ffd::chebyshev_polynomial_s::matsubara_transform<n_t>(0);
  ffd::l_array<Complex, size_x_y> YRZ_si2_ko0, BCS_si2_ko0, BCD_si2_ko0, mix_si2_ko0,
    YRZ_si3_ko0, BCS_si3_ko0, BCD_si3_ko0, mix_si3_ko0, FLAT_si2_ko0, FLAT_si3_ko0;
  for(std::size_t k=0; k<size_x_y; ++k){
    YRZ_si2_ko0[k] = 0.;
    for(std::size_t j=0; j<n_t; ++j){
      YRZ_si2_ko0[k] += Beta*YRZ_si2_kt[k].coef[j]*ch_to_mts[j];
    }


    YRZ_si3_ko0[k] = 0.;
    for(std::size_t j=0; j<n_t; ++j){
      YRZ_si3_ko0[k] += Beta*YRZ_si3_kt[k].coef[j]*ch_to_mts[j];
    }


    BCS_si2_ko0[k] = 0.;
    for(std::size_t j=0; j<n_t; ++j){
      BCS_si2_ko0[k] += Beta*BCS_si2_kt[k].coef[j]*ch_to_mts[j];
    }


    BCS_si3_ko0[k] = 0.;
    for(std::size_t j=0; j<n_t; ++j){
      BCS_si3_ko0[k] += Beta*BCS_si3_kt[k].coef[j]*ch_to_mts[j];
    }


    BCD_si2_ko0[k] = 0.;
    for(std::size_t j=0; j<n_t; ++j){
      BCD_si2_ko0[k] += Beta*BCD_si2_kt[k].coef[j]*ch_to_mts[j];
    }

    
    mix_si2_ko0[k] = 0.;
    for(std::size_t j=0; j<n_t; ++j){
      mix_si2_ko0[k] += Beta*mix_si2_kt[k].coef[j]*ch_to_mts[j];
    }


    mix_si3_ko0[k] = 0.;
    for(std::size_t j=0; j<n_t; ++j){
      mix_si3_ko0[k] += Beta*mix_si3_kt[k].coef[j]*ch_to_mts[j];
    }

    
    FLAT_si2_ko0[k] = 0.;
    for(std::size_t j=0; j<n_t; ++j){
      FLAT_si2_ko0[k] += Beta*FLAT_si2_kt[k].coef[j]*ch_to_mts[j];
    }


    FLAT_si3_ko0[k] = 0.;
    for(std::size_t j=0; j<n_t; ++j){
      FLAT_si3_ko0[k] += Beta*FLAT_si3_kt[k].coef[j]*ch_to_mts[j];
    }

  }


  std::ofstream file_spc_0("out_data_temp/spc_0.dat");
  std::ofstream file_G_0("out_data_temp/G_0.dat");
  
  std::ofstream file_spc_YRZ("out_data_temp/spc_YRZ.dat");
  std::ofstream file_isi_YRZ("out_data_temp/isi_YRZ.dat");
  std::ofstream file_rsi_YRZ("out_data_temp/rsi_YRZ.dat");
  std::ofstream file_spc2_YRZ("out_data_temp/spc2_YRZ.dat");
  std::ofstream file_isi2_YRZ("out_data_temp/isi2_YRZ.dat");
  std::ofstream file_rsi2_YRZ("out_data_temp/rsi2_YRZ.dat");
  std::ofstream file_spc3_YRZ("out_data_temp/spc3_YRZ.dat");
  std::ofstream file_isi3_YRZ("out_data_temp/isi3_YRZ.dat");
  std::ofstream file_rsi3_YRZ("out_data_temp/rsi3_YRZ.dat");
  std::ofstream file_G_YRZ("out_data_temp/G_YRZ.dat");
  
  std::ofstream file_spc_BCS("out_data_temp/spc_BCS.dat");
  std::ofstream file_isi_BCS("out_data_temp/isi_BCS.dat");
  std::ofstream file_rsi_BCS("out_data_temp/rsi_BCS.dat");
  std::ofstream file_spc2_BCS("out_data_temp/spc2_BCS.dat");
  std::ofstream file_isi2_BCS("out_data_temp/isi2_BCS.dat");
  std::ofstream file_rsi2_BCS("out_data_temp/rsi2_BCS.dat");
  std::ofstream file_spc3_BCS("out_data_temp/spc3_BCS.dat");
  std::ofstream file_isi3_BCS("out_data_temp/isi3_BCS.dat");
  std::ofstream file_rsi3_BCS("out_data_temp/rsi3_BCS.dat");
  

  std::ofstream file_G_BCS("out_data_temp/G_BCS.dat");
  
  std::ofstream file_spc_BCD("out_data_temp/spc_BCD.dat");
  std::ofstream file_isi_BCD("out_data_temp/isi_BCD.dat");
  std::ofstream file_rsi_BCD("out_data_temp/rsi_BCD.dat");
  std::ofstream file_spc2_BCD("out_data_temp/spc2_BCD.dat");
  std::ofstream file_isi2_BCD("out_data_temp/isi2_BCD.dat");
  std::ofstream file_rsi2_BCD("out_data_temp/rsi2_BCD.dat");
  std::ofstream file_G_BCD("out_data_temp/G_BCD.dat");

  std::ofstream file_spc_mix("out_data_temp/spc_mix.dat");
  std::ofstream file_isi_mix("out_data_temp/isi_mix.dat");
  std::ofstream file_rsi_mix("out_data_temp/rsi_mix.dat");
  std::ofstream file_spc2_mix("out_data_temp/spc2_mix.dat");
  std::ofstream file_isi2_mix("out_data_temp/isi2_mix.dat");
  std::ofstream file_rsi2_mix("out_data_temp/rsi2_mix.dat");
  std::ofstream file_G_mix("out_data_temp/G_mix.dat");
  
  std::ofstream file_spc_FLAT("out_data_temp/spc_FLAT.dat");
  std::ofstream file_isi_FLAT("out_data_temp/isi_FLAT.dat");
  std::ofstream file_rsi_FLAT("out_data_temp/rsi_FLAT.dat");
  std::ofstream file_spc2_FLAT("out_data_temp/spc2_FLAT.dat");
  std::ofstream file_isi2_FLAT("out_data_temp/isi2_FLAT.dat");
  std::ofstream file_rsi2_FLAT("out_data_temp/rsi2_FLAT.dat");
  std::ofstream file_G_FLAT("out_data_temp/G_FLAT.dat");
    std::ofstream file_isi3_FLAT("out_data_temp/isi3_FLAT.dat");
  std::ofstream file_rsi3_FLAT("out_data_temp/rsi3_FLAT.dat");

  
  Real G0_max = 0., BCS_G_max = 0., PARAB_max=0.;
  for(std::size_t nx=0; nx<=size_x/2; ++nx){
    for(std::size_t ny=0; ny<=size_y/2; ++ny){
      Real const kx = 2.*nx*Pi/size_x;
      Real const ky = 2.*ny*Pi/size_y;

      file_spc_0<<kx<<" "<<ky<<" "<<
	-imag(1./(I*omg_mts(0) - xi_k(kx, ky) + mu_0 ))/Pi<<'\n';
      file_G_0<<nx<<" "<<ny<<" "<<
	-G0_rt[ny+size_y*nx](Beta)<<'\n';
      if(auto g0_max = std::abs((1./(I*omg_mts(0) - xi_k(kx, ky) + mu_0 )));
	 g0_max > G0_max){
	G0_max = g0_max;

      }
      
      Complex const si_YRZ = YRZsi_ko(kx, ky, 0);
      file_isi_YRZ<<kx<<" "<<ky<<" "<<imag(si_YRZ)<<"\n";
      file_rsi_YRZ<<kx<<" "<<ky<<" "<<real(si_YRZ)<<"\n";
      file_spc_YRZ<<kx<<" "<<ky<<" "<<
	-imag(1./(I*omg_mts(0) - xi_k(kx, ky) + mu_0 - si_YRZ))/Pi<<'\n';
      Complex const si2_YRZ = YRZ_si2_ko0[ny+size_y*nx];
      file_isi2_YRZ<<kx<<" "<<ky<<" "<<imag(si2_YRZ)<<"\n";
      file_rsi2_YRZ<<kx<<" "<<ky<<" "<<real(si2_YRZ)<<"\n";
      file_spc2_YRZ<<kx<<" "<<ky<<" "<<
	-imag(1./(I*omg_mts(0) - xi_k(kx, ky) + mu_0 - si2_YRZ))/Pi<<'\n';
      file_G_YRZ<<nx<<" "<<ny<<" "<<
	-YRZ_G_rt[ny+size_y*nx](Beta)<<'\n';


      Complex const si_BCS = BCSsi_ko_mu0(kx, ky, 0, mu_0);
      file_isi_BCS<<kx<<" "<<ky<<" "<<imag(si_BCS)<<"\n";
      file_rsi_BCS<<kx<<" "<<ky<<" "<<real(si_BCS)<<"\n";
      file_spc_BCS<<kx<<" "<<ky<<" "<<
	-imag(1./(I*omg_mts(0) - xi_k(kx, ky) + mu_0 - si_BCS))/Pi<<'\n';
      Complex const si2_BCS = BCS_si2_ko0[ny+size_y*nx];
      file_isi2_BCS<<kx<<" "<<ky<<" "<<imag(si2_BCS)<<"\n";
      file_rsi2_BCS<<kx<<" "<<ky<<" "<<real(si2_BCS)<<"\n";
      file_spc2_BCS<<kx<<" "<<ky<<" "<<
	-imag(1./(I*omg_mts(0) - xi_k(kx, ky) + mu_0 - si2_BCS))/Pi<<'\n';
      file_G_BCS<<nx<<" "<<ny<<" "<<
	-BCS_G_rt[ny+size_y*nx](Beta)<<'\n';
      if(auto bcs_max=std::abs((1./(I*omg_mts(0) - xi_k(kx, ky) + mu_0 - si_BCS)));
	 bcs_max>BCS_G_max){
	BCS_G_max = bcs_max;
      }


      Complex const si_BCD = BCDsi_ko_mu0(kx, ky, 0, mu_0);
      file_isi_BCD<<kx<<" "<<ky<<" "<<imag(si_BCD)<<"\n";
      file_rsi_BCD<<kx<<" "<<ky<<" "<<real(si_BCD)<<"\n";
      file_spc_BCD<<kx<<" "<<ky<<" "<<
	-imag(1./(I*omg_mts(0) - xi_k(kx, ky) + mu_0 - si_BCD))/Pi<<'\n';
      Complex const si2_BCD = BCD_si2_ko0[ny+size_y*nx];
      file_isi2_BCD<<kx<<" "<<ky<<" "<<imag(si2_BCD)<<"\n";
      file_rsi2_BCD<<kx<<" "<<ky<<" "<<real(si2_BCD)<<"\n";
      file_spc2_BCD<<kx<<" "<<ky<<" "<<
	-imag(1./(I*omg_mts(0) - xi_k(kx, ky) + mu_0 - si2_BCD))/Pi<<'\n';
      file_G_BCD<<nx<<" "<<ny<<" "<<
	-BCD_G_rt[ny+size_y*nx](Beta)<<'\n';


      Complex const si_mix = PARABsi_ko_mu0(kx, ky, 0, mu_0);
	// mix_BCS*BCSsi_ko_mu0(kx, ky, 0, mu_0)+
	//     mix_BCD*BCDsi_ko_mu0(kx, ky, 0, mu_0)+mix_YRZ*YRZsi_ko(kx, ky, 0);;
      file_isi_mix<<kx<<" "<<ky<<" "<<imag(si_mix)<<"\n";
      file_rsi_mix<<kx<<" "<<ky<<" "<<real(si_mix)<<"\n";
      file_spc_mix<<kx<<" "<<ky<<" "<<
	-imag(1./(I*omg_mts(0) - xi_k(kx, ky) + mu_0 - si_mix))/Pi<<'\n';
      Complex const si2_mix = mix_si2_ko0[ny+size_y*nx];
      file_isi2_mix<<kx<<" "<<ky<<" "<<imag(si2_mix)<<"\n";
      file_rsi2_mix<<kx<<" "<<ky<<" "<<real(si2_mix)<<"\n";
      file_spc2_mix<<kx<<" "<<ky<<" "<<
	-imag(1./(I*omg_mts(0) - xi_k(kx, ky) + mu_0 - si2_mix))/Pi<<'\n';
      file_G_mix<<nx<<" "<<ny<<" "<<
	-mix_G_rt[ny+size_y*nx](Beta)<<'\n';
      if(auto bcs_max=std::abs((1./(I*omg_mts(0) - xi_k(kx, ky) + mu_0 - si_mix)));
	 bcs_max>PARAB_max){
	PARAB_max = bcs_max;
      }


      Complex const si_FLAT = PARABsi_ko_mu0(kx, ky, 0, mu_0);
	// FLAT_BCS*BCSsi_ko_mu0(kx, ky, 0, mu_0)+
	//     FLAT_BCD*BCDsi_ko_mu0(kx, ky, 0, mu_0)+FLAT_YRZ*YRZsi_ko(kx, ky, 0);;
      file_isi_FLAT<<kx<<" "<<ky<<" "<<imag(si_FLAT)<<"\n";
      file_rsi_FLAT<<kx<<" "<<ky<<" "<<real(si_FLAT)<<"\n";
      file_spc_FLAT<<kx<<" "<<ky<<" "<<
	-imag(1./(I*omg_mts(0) - xi_k(kx, ky) + mu_0 - si_FLAT))/Pi<<'\n';
      Complex const si2_FLAT = FLAT_si2_ko0[ny+size_y*nx];
      file_isi2_FLAT<<kx<<" "<<ky<<" "<<imag(si2_FLAT)<<"\n";
      file_rsi2_FLAT<<kx<<" "<<ky<<" "<<real(si2_FLAT)<<"\n";
      file_spc2_FLAT<<kx<<" "<<ky<<" "<<
	-imag(1./(I*omg_mts(0) - xi_k(kx, ky) + mu_0 - si2_FLAT))/Pi<<'\n';
      file_G_FLAT<<nx<<" "<<ny<<" "<<
	-FLAT_G_rt[ny+size_y*nx](Beta)<<'\n';
      if(auto bcs_max=std::abs((1./(I*omg_mts(0) - xi_k(kx, ky) + mu_0 - si_FLAT)));
	 bcs_max>PARAB_max){
	PARAB_max = bcs_max;
      }

      Complex const si3_FLAT = FLAT_si3_ko0[ny+size_y*nx];
      file_isi3_FLAT<<kx<<" "<<ky<<" "<<imag(si3_FLAT)<<"\n";
      file_rsi3_FLAT<<kx<<" "<<ky<<" "<<real(si3_FLAT)<<"\n";
  
    }
  }
  std::cerr<<std::setprecision(10)<<"G0_max = "<<G0_max<<std::endl;

	std::cerr<<"BCS_G_max = "<<BCS_G_max<<std::endl;
	std::cerr<<"PARAB_G_max = "<<PARAB_max<<std::endl;
  { 
    auto Lx = (1ul<<l_x);
    auto Ly = (1ul<<l_y);
    std::ostringstream file_sigma_name, file_G_name, common_part,
      run_name, common_part_2;
    common_part_2<<"out_data/";
    run_name<<"BCS_mu"<<mu<<"_Dlt"<<BCSdlt;
    if(std::abs(BCSdlt) > 1e-15){
      run_name<<"_B"<<BCSBeta;
    }
    common_part<<std::setprecision(10)<<"_B"<<Beta;
    common_part<<"_U"<<U;
    common_part<<"_"<<Lx<<"x"<<Ly;
    common_part<<"_Nc"<<N_Chebyshev_input;
    common_part<<"_AS"<<alpha_shift;
    common_part<<"_name"<<run_name.str()<<".dat";
    file_sigma_name<<common_part_2.str()<<"SigmaR";
    file_G_name    <<common_part_2.str()<<"GR____";
    file_sigma_name<<common_part.str();
    file_G_name<<common_part.str();


    std::ofstream fileGBCS(file_G_name.str());
    std::ofstream fileSBCS(file_sigma_name.str());



    auto const chebyshev_nodes = ffd::chebyshev_polynomial_s::
    create_chebyshev_nodes<n_t>({0., Beta});
    for( std::size_t x=0; x < (1<<l_x); ++x){
      for( std::size_t y=0; y < (1<<l_y); ++y){
	std::size_t sp_ind = x*(1<<l_y)+y;
	for( std::size_t k=0; k < n_t; ++k){
	  fileGBCS<<std::setprecision(12)<<BCS_G_rt[sp_ind](chebyshev_nodes[k])<<'\n';
	  fileSBCS<<std::setprecision(12)<<BCS_si_rt[sp_ind](chebyshev_nodes[k])<<'\n';
	}
      }
    }
  }




   {
    auto Lx = (1ul<<l_x);
    auto Ly = (1ul<<l_y);
    std::ostringstream file_sigma_name, file_G_name, common_part,
      run_name, common_part_2;
    common_part_2<<"out_data/";
    run_name<<"mix_mu"<<mu<<"_Dlt"<<PARABdlt;
    common_part<<std::setprecision(10)<<"_B"<<Beta;
    common_part<<"_U"<<U;
    common_part<<"_"<<Lx<<"x"<<Ly;
    common_part<<"_Nc"<<N_Chebyshev_input;
    common_part<<"_AS"<<alpha_shift;
    common_part<<"_name"<<run_name.str()<<".dat";
    file_sigma_name<<common_part_2.str()<<"SigmaR";
    file_G_name    <<common_part_2.str()<<"GR____";
    file_sigma_name<<common_part.str();
    file_G_name<<common_part.str();


    std::ofstream fileGPARAB(file_G_name.str());
    std::ofstream fileSPARAB(file_sigma_name.str());



    auto const chebyshev_nodes = ffd::chebyshev_polynomial_s::
    create_chebyshev_nodes<n_t>({0., Beta});
    for( std::size_t x=0; x < (1<<l_x); ++x){
      for( std::size_t y=0; y < (1<<l_y); ++y){
	std::size_t sp_ind = x*(1<<l_y)+y;
	for( std::size_t k=0; k < n_t; ++k){
	  fileGPARAB<<std::setprecision(12)<<mix_G_rt[sp_ind](chebyshev_nodes[k])<<'\n';
	  fileSPARAB<<std::setprecision(12)<<mix_si_rt[sp_ind](chebyshev_nodes[k])<<'\n';
	}
      }
    }
  }



   
   {
    auto Lx = (1ul<<l_x);
    auto Ly = (1ul<<l_y);
    std::ostringstream file_sigma_name, file_G_name, common_part,
      run_name, common_part_2;
    common_part_2<<"out_data/";
    run_name<<"FLAT_mu"<<mu<<"_tR"<<tR<<"_dmuR"<<dmuR;
    common_part<<std::setprecision(10)<<"_B"<<Beta;
    common_part<<"_U"<<U;
    common_part<<"_"<<Lx<<"x"<<Ly;
    common_part<<"_Nc"<<N_Chebyshev_input;
    common_part<<"_AS"<<alpha_shift;
    common_part<<"_name"<<run_name.str()<<".dat";
    file_sigma_name<<common_part_2.str()<<"SigmaR";
    file_G_name    <<common_part_2.str()<<"GR____";
    file_sigma_name<<common_part.str();
    file_G_name<<common_part.str();


    std::ofstream fileGFLAT(file_G_name.str());
    std::ofstream fileSFLAT(file_sigma_name.str());



    auto const chebyshev_nodes = ffd::chebyshev_polynomial_s::
    create_chebyshev_nodes<n_t>({0., Beta});
    for( std::size_t x=0; x < (1<<l_x); ++x){
      for( std::size_t y=0; y < (1<<l_y); ++y){
	std::size_t sp_ind = x*(1<<l_y)+y;
	for( std::size_t k=0; k < n_t; ++k){
	  fileGFLAT<<std::setprecision(12)<<FLAT_G_rt[sp_ind](chebyshev_nodes[k])<<'\n';
	  fileSFLAT<<std::setprecision(12)<<FLAT_si_rt[sp_ind](chebyshev_nodes[k])<<'\n';
	}
      }
    }
  }

}
