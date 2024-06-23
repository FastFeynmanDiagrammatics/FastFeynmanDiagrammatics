namespace ffd::gauss_kronrod{  
  struct ConstexprIfFailure{template<typename> void failure();};
  

  template<int j,
	   typename func_t>
  
  std::array<Real, 2>
  
  Integrate_s(func_t&& func_v,
	      std::array<Real, 2> interval){
    Real const mean = .5*( interval[1] + interval[0] );
    Real const diff = .5*( interval[1] - interval[0] );


    Real gauss = 0., kronrod = 0.;
    for( std::size_t u = 0; u < j; ++u ){
      if constexpr( j == 3 ){
	  Real const f_eval = func_v(mean + diff*N3[u]);
	  kronrod += f_eval*K3[u];
	  auto const& is_odd = u&1;
	  gauss += is_odd*f_eval*G3[is_odd*(u>>1)];
	}else if constexpr( j == 7 ){
	  Real const f_eval = func_v(mean + diff*N7[u]);
	  kronrod += f_eval*K7[u];
	  auto const& is_odd = u&1;
	  gauss += is_odd*f_eval*G7[is_odd*(u>>1)];
	}else if constexpr( j == 15 ){
	  Real const f_eval = func_v(mean + diff*N15[u]);
	  kronrod += f_eval*K15[u];
	  auto const& is_odd = u&1;
	  gauss += is_odd*f_eval*G15[is_odd*(u>>1)];
	}else if constexpr( j == 31 ){
	  Real const f_eval = func_v(mean + diff*N31[u]);
	  kronrod += f_eval*K31[u];
	  auto const& is_odd = u&1;
	  gauss += is_odd*f_eval*G31[is_odd*(u>>1)];
	}else if constexpr( j == 63 ){
	  Real const f_eval = func_v(mean + diff*N63[u]);
	  kronrod += f_eval*K63[u];
	  auto const& is_odd = u&1;
	  gauss += is_odd*f_eval*G63[is_odd*(u>>1)];
	}else{
	ConstexprIfFailure failure; failure.failure<bool>();
      }
    }
    
    
    std::array<Real, 2> ret{kronrod*diff, std::abs((kronrod-gauss)*diff)};
    return ret;
  }


  
  template<typename func_t>
  
  Real
  
  Integrate_63(func_t&& func_v,
	       std::array<Real, 2> interval){
    Real const mean = .5*( interval[1] + interval[0] );
    Real const diff = .5*( interval[1] - interval[0] );


    Real kronrod = 0.;
    for(std::size_t u=0; u<63; ++u){
      kronrod += K63[u]*func_v(mean + diff*N63[u]);;
    }
    
    
    return kronrod*diff;
  }


}//namespace
