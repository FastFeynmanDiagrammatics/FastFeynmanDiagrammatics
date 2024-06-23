namespace ffd::user_space{

  Real G0(Real tau){
    using std::exp;
    Real sign = 1.;
    while( tau < 0 ){tau += Beta; sign = -sign; }
    while( tau > Beta ){ tau -= Beta; sign = -sign; }
    return -sign*exp(tau*mu)/(1+exp(Beta*mu));
  }


  Real matsu_freq(int omega_n){
    using ffd::core_math::Pi;
    return (2*omega_n+1)*Pi/Beta;
  }


  Complex G0_omega(int omega_n){
    return 1./(I*matsu_freq(omega_n)+mu);
  }


  Complex Sigma_R_omega(int omega_n){
    return .25*Delta*Delta*U*U/(I*matsu_freq(omega_n)-mu*Beta_Sigma/Beta);
  }

  
  Complex GR_omega(int omega_n){
    return 1./(1./G0_omega(omega_n)-Sigma_R_omega(omega_n));
  }


  Real GR(Real tau){
    using std::exp;
    Complex delta_GR = 0.;
    for(int n_freq = -n_frequencies_DMFT; n_freq<0; n_freq++){
      delta_GR += (GR_omega(n_freq)-G0_omega(n_freq))*exp(-I*matsu_freq(n_freq)*tau);
    }
    delta_GR /= Beta;
    return 2*std::real(delta_GR)+G0(tau);
  }
}
