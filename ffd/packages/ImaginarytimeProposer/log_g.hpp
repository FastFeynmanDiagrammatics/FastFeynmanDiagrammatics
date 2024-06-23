namespace ffd::itime_proposer{

  auto
  log_g(Real Beta,
	Real t0){

   using ffd::user_space::Proba;	
	  
    Real const A =
      .5/(t0*std::log(1+.5*Beta/t0));

    Real const A_t =
      1./(A*t0);

    auto Periodize =
      ffd::periodic_numbers::Periodic_g({0., Beta});

    auto Periodize_2 =
      ffd::periodic_numbers::Periodic_g({-.5*Beta, .5*Beta});

    auto propose =
      [A_t, t0, Periodize]
      (Real t){
	Real const p = .5*Proba();
	Real const tp = t0*(std::exp(p*A_t)-1.);
	return Periodize(t+(1-2*(Proba()>.5))*tp).first;
      };
    
    auto conditional =
      [A, t0, Periodize_2]
      (Real tau0, Real tau1){
	auto const [t01, discard] = Periodize_2(tau0-tau1);
	return A/(1+std::abs(t01)/t0);
      };


    auto ret =
      [conditional, propose, Beta]
      (auto&& ...x){
	if constexpr(sizeof...(x) == 0){
	    return .5*Beta;
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
