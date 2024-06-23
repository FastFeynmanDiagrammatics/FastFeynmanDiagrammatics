  auto d_wave_gap_k =
    [=](Real kx, Real ky){
      return cos(kx)-cos(ky);
    };
  auto xi_HF_k =
    [=](Real kx, Real ky){
      return -2.*t_hopping*(cos(kx)+cos(ky));
    };
  auto xi_k =
    [=](Real kx, Real ky){
      Real const cos_kx = cos(kx);
      Real const cos_ky = cos(ky);
      return -2.*t_hopping*(cos_kx+cos_ky) -
	4.*tp_hopping*cos_kx*cos_ky;
    };
  auto omg_mts =
    [=](int omg){
      return Pi*(1.+2.*omg)/Beta;
    };
  auto green_function_Et =
    [=](Real E, Real t){
      Real constexpr Beta_2 = .5*Beta;
      Real const EB_2 = E*Beta_2;
      return -.5*exp(EB_2-E*t)/cosh(EB_2);
    };
  auto green_function_Eo =
    [=](Real E, int omg){
      return 1./(I*omg_mts(omg)-E);
    };
  auto YRZsi_kt =
    [=](Real kx, Real ky, Real t){
      Real YRZ2 = YRZdlt*d_wave_gap_k(kx, ky);
      YRZ2 *= YRZ2;
      return YRZ2*green_function_Et(-YRZBeta*xi_HF_k(kx, ky)/Beta, t);
    };
  auto BCDsi_kt_mu0 =
    [=](Real kx, Real ky, Real t, Real mu0){
      Real BCD2 = BCDdlt*d_wave_gap_k(kx, ky);
      BCD2 *= BCD2;
      return BCD2*green_function_Et(-BCDBeta*(xi_k(kx, ky)-mu0)/Beta, t);
    };
  auto BCSsi_kt_mu0 =
    [=](Real kx, Real ky, Real t, Real mu0){
      return BCSdlt*BCSdlt*green_function_Et(-BCSBeta*(xi_k(kx, ky)-mu0)/Beta-a_sh_BCS, t);
    };
  auto PARABsi_kt_mu0 =
    [=](Real kx, Real ky, Real t, Real mu0){
      Real x = PARABdlt;
      x *= x;
      x *= x;
      return -.25*x*t*(Beta-t);
    };
  auto FLATsi_kt_mu0 =
    [=](Real kx, Real ky, Real t, Real mu0){
      return -.5*(-(1.-tR)*(xi_k(kx, ky)-mu0)-dmuR)*pow(omg_mts(0), 2)*(t-.5*Beta);
    };
  // auto FLATsi_kt_mu0 =
  //   [=](Real kx, Real ky, Real t, Real mu0){
  //     return -.5*((1.-tR)*xi_k(kx, ky)-dmuR)*pow(omg_mts(0), 2)*(t-.5*Beta);
  //   };

  auto green_function0_kt_mu0 =
    [=](Real kx, Real ky, Real t, Real mu0){
      return green_function_Et(xi_k(kx, ky)-mu0, t);
    };
  auto green_function0_ko_mu0 =
    [=](Real kx, Real ky, int omg, Real mu0){
      return green_function_Eo(xi_k(kx, ky)-mu0, omg);
    };
  auto YRZsi_ko =
    [=](Real kx, Real ky, int omg){
      Real YRZ2 = YRZdlt*d_wave_gap_k(kx, ky);
      YRZ2 *= YRZ2;
      return YRZ2*green_function_Eo(-YRZBeta*xi_HF_k(kx, ky)/Beta, omg);
    };
  auto BCSsi_ko_mu0 =
    [=](Real kx, Real ky, int omg, Real mu0){
      return BCSdlt*BCSdlt*green_function_Eo(-BCSBeta*(xi_k(kx, ky)-mu0)/Beta-a_sh_BCS, omg);
    };
  auto BCDsi_ko_mu0 =
    [=](Real kx, Real ky, int omg, Real mu0){
      Real BCD2 = BCDdlt*d_wave_gap_k(kx, ky);
      BCD2 *= BCD2;
      return BCD2*green_function_Eo(-BCDBeta*(xi_k(kx, ky)-mu0)/Beta, omg);
    };
  auto PARABsi_ko_mu0 = 
    [=](Real kx, Real ky, int omg, Real mu0){
      Real x = PARABdlt;
      x *= x;
      x *= x;
      return x/(I*std::pow(omg_mts(omg), 3));
    };
  auto FLATsi_ko_mu0 =
    [=](Real kx, Real ky, int omg, Real mu0){
      return (-(1.-tR)*(xi_k(kx, ky)-mu0)-dmuR)*pow(omg_mts(0)/omg_mts(omg), 2);
    };
   // auto FLAT2si_omg =
   //  [=](int omg){
   //    Real const omg_0 = omg_mts(0);
   //    Real const omg_n = omg_mts(omg);
   //    Real x_n = omg_n/LambdaR_omg_flat;
   //    x_n *= x_n;
   //    Real x_0 = omg_0/LambdaR_omg_flat;
   //    x_0 *= x_0;
   //    return exp(-.5*(x_n-x_0));
   //  };
  auto FLAT2si_omg =
    [=](int omg){
      Real const omg_0 = omg_mts(0);
      Real const omg_n = omg_mts(omg);
      return std::pow(omg_0/omg_n, 2);
    };
 
  auto FLAT2si_t =
    [=](Real t){
      Real re = 0;
      for(int j=-n_freq_mts; j<=n_freq_mts; ++j){
	Real const omega_t = omg_mts(j)*t;
	re += FLAT2si_omg(j)*cos(omega_t);
      }
      return re/Beta;
    };
  auto FLATsi_k_mu0 =
    [=](Real kx, Real ky, Real mu0){
      return -(1.-tR)*(xi_k(kx, ky)-mu0) - dmuR;
    };
  auto FLAT2si_k_mu0 =
    [=](Real kx, Real ky, Real mu0){
      Real xik = xi_k(kx, ky)-mu0;
      return FLATsi_k_mu0(kx, ky, mu0)/(1+0.*xik*xik);
    };

  ffd::chebyshev_polynomial_s::TPoly<n_t> flat2_ch_t(FLAT2si_t, {0., Beta});
auto FLAT2si_kt_mu0 =
    [=](Real kx, Real ky, Real t, Real mu0){
      return FLAT2si_k_mu0(kx, ky, mu0)*flat2_ch_t(t);
    };
