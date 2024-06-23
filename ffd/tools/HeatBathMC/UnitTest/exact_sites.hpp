namespace ffd::heat_bath_mc::unit_test{

  template<int order,
	   int order_plus = 1,
	   int half_num_sites = 1,
	   bool sum_over_config = true>
  void
  exact_sites(){

    auto constexpr num_sites = 2*half_num_sites+1;
    auto N_iter = 1ul<<20;

    
    auto f = [](auto r){
	       Real arg_exp = 0.;
	       int sign = 1;
	       for(int j=0; j<size(r); ++j){
		 arg_exp += r[j]*r[j];
		 if( (r[j]&1) == 1){
		   sign = - sign;
		 }
	       }
	       return sign/(1.+2*arg_exp);
		 //std::exp(-std::sqrt(arg_exp))*sign;
	     };

    
    auto phi_unnorm = [&](int x){
			x = x<0 ? -x : x;
			if(x>half_num_sites){
			  x -= num_sites;
			}
			return 1./(1.+std::abs(x));
		      };

    
    Real norm_phi = 0.;
    for(int j=0; j<num_sites; ++j){
      norm_phi += phi_unnorm(j);
    }


    auto phi = [&](int x){
		 return phi_unnorm(x)/norm_phi;
	       };
    

    auto pseudo_pot = [&](auto r, auto r_added){
			Real pseudo_potential = 0.;
			for( std::size_t j=0; j < order; ++j){
			  pseudo_potential += phi(r[j]-r_added);
			}
			return pseudo_potential/order;
		      };


    auto proposer =
      [&](int r){
	Real proba = ffd::user_space::Proba();
	Real proba_cumul = 0.;
	int j;
	for(j=0; j<num_sites; ++j){
	  proba_cumul += phi(j);
	  if(proba_cumul>proba){
	    break;
	  }
	}
	return (r+j)%num_sites;
      };
    
    
    std::array<int, order> r;
    std::array<int, order+order_plus> R;
    BinaryInt const cardinality_order_plus = (1ul<<(order+order_plus));
    ffd::accumulator_markov_chain::
      AccumulatorMarkovChain average_sign;
    r.fill(0);
    for( long t_MC = 0; t_MC < N_iter; ++t_MC){
      for(int j=0; j<order; ++j){
	R[j] = r[j];
      }
      for(int j=0; j<order_plus;++j){
	int pos = ffd::random_distributions::RandomInRange(order);
	R[order+j] = proposer(r[pos]);
      }


      std::vector<Real> fluxes;
      std::vector<BinaryInt> fluxes_S;
      std::vector<Real> sign_S;
      for( BinaryInt S=0; S < cardinality_order_plus; ++S ){
	if( __builtin_popcount(S) == order ){
	  auto vec_S = ffd::set_theory::VectorOfBinaryDigitsOf(S);
	  auto vec_compl_S = ffd::set_theory::
	    VectorOfBinaryDigitsOf(cardinality_order_plus - 1 - S);
	  std::array<int, order> r_temp;
	  std::array<int, order_plus> r_compl_temp;
	  for( int j=0; j < order; ++j){
	    r_temp[j] = R[vec_S[j]];
	  }
	  for( int j=0; j < order_plus; ++j){
	    r_compl_temp[j] = R[vec_compl_S[j]];
	  }

	  
	  fluxes_S.push_back(S);
	  auto const f_val = f(r_temp);
	  Real pseudo_potential = 1.;
	  for(int j=0; j<order_plus; ++j){
	    pseudo_potential *= pseudo_pot(r_temp, r_compl_temp[j]);
	  }
	  fluxes.push_back(pseudo_potential*f_val);
	  sign_S.push_back(f_val/std::abs(f_val));
	}
      }


      auto [time_S, discard, pos_conf] =
	heat_bath_step(fluxes, std::vector<Real>());

      
      auto const vec_S = ffd::set_theory::
	VectorOfBinaryDigitsOf(fluxes_S[pos_conf]);
      for( int j=0; j<order; ++j){
	r[j] = R[vec_S[j]];
      }
      
      if constexpr(sum_over_config){
	  Real temp_sign = 0.;
	    for( std::size_t j=0; j<size(fluxes); ++j){
	      temp_sign +=
		time_S[j]*sign_S[j];
	    }
	    average_sign += temp_sign;
	}else{
	average_sign += sign_S[pos_conf];
      }
    }


    auto [mean, error] = average_sign.MeanValueError();
    std::cerr<<"order = "<<order<<", order_plus = "<<order_plus<<", L = "<<num_sites<<std::endl;
    std::cerr<<"<f/|f|> = "<<mean<<"+/-"<<error<<", \u03C4 = "<<average_sign.AutocorrelationTime()<<", var = "<<average_sign.Variance()<<std::endl;


    std::array<int, order> L;
    L.fill(num_sites);
    Real exact_f = 0.;
    Real exact_af = 0.;
    for( auto r: ffd::user_space::CartesianProductRange(L) ){
      auto f_ex = f(r);
      exact_f += f_ex;
      exact_af += std::abs(f_ex);
    }


    std::cerr<<"<f/|f|> = "<<exact_f/exact_af<<std::endl;
  }

}//namespace
