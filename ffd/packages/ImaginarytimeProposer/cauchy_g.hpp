namespace ffd::itime_proposer{

  // zero temperature version
  auto
  cauchy_g(Real t0){

   using ffd::user_space::Proba;
   using ffd::core_math::Pi;
	  
    auto propose =
      [t0]
      (Real t){
	Real const p = Proba();
	Real const tp = t0*std::tan(Pi*p);
	if(Proba()<.5){
	  return t + tp;
	}else{
	  return t - tp;
	}
      };

    Real const A = 1./(Pi*t0);
    
    auto conditional =
      [t0, A]
      (Real tau0, Real tau1){
	Real const dtau = std::abs(tau0-tau1)/t0;
	return A/(1+dtau*dtau);
      };

    
    auto ret =
      [conditional, propose, t0]
      (auto&& ...x){
	if constexpr(sizeof...(x) == 0){
	    return 0.;
	  }else{
	  if constexpr(sizeof...(x) == 1){
	      return propose(std::forward<decltype(x)>(x)...);
	    }else{
	    return conditional(std::forward<decltype(x)>(x)...);
	  }
	}
      };


    return ret;
  }

}//namespace
