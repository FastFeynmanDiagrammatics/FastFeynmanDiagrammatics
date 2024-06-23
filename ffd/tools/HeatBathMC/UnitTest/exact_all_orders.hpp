namespace ffd::heat_bath_mc::unit_test{

  template<int order,
	   int order_plus,
	   int half_num_sites,
	   bool is_translation_invariant,
	   typename lambda_function_t>
  
  void exact_all_orders(lambda_function_t f,
			unsigned long N_iter = 1ul<<18){

    auto constexpr num_sites = 2*half_num_sites+1;

    
    auto phi_unnorm = [=](int x){
			x = x>=0 ? x : -x;
			if(x>half_num_sites){
			  x -= num_sites;
			}
			return 1./(1.+x*x*x*x);
		      };

    
    Real norm_phi = 0.;
    for(int j=0; j<num_sites; ++j){
      norm_phi += phi_unnorm(j);
    }


    auto phi = [=](int x1, int x2){
		 return phi_unnorm(x1-x2)/norm_phi;
	       };

    auto phi_proposer = [=](int x){
		      Real p = ffd::user_space::Proba();
		      int j=0;
		      Real p_cumul = 0.;
		      for(;; ++j){
			p_cumul += phi(x, j);
			if(p_cumul > p){
			  break;
			}
		      }
		      return j;
		    };

    
    auto Phi = [=](auto... x){
		 if constexpr(sizeof...(x) == 2){
		     return phi(x...);
		   }else if constexpr(sizeof...(x) == 1){
		     return phi_proposer(x...);
		   }else{
		   assert(false);
		   return 0;
		 }
	       };


    auto compute_P =
      [=](auto x){
	std::array<Real, (1<<(order+order_plus))> P;
	P.fill(0.);
	P[0] = 1.;
	for(BinaryInt S=1; S<(1<<(order+order_plus)); ++S){
	  // if(__builtin_popcount(S) <= order){
	  std::vector<BinaryInt> const digits_S = ffd::set_theory::VectorOfBinaryDigitsOf(S);
	  std::size_t const card_S = size(digits_S);
	  std::vector<int> x_reduced(card_S);
	  for(std::size_t j=0; j<card_S; ++j){
	    x_reduced[j] = x[digits_S[j]];
	  }
	  P[S] = f(x_reduced);
	// }
	}
	return P;
      };

    
    std::array<int, order+order_plus> X;
    std::array<ffd::accumulator_markov_chain::
	       AccumulatorMarkovChain<>, order+order_plus+1> orders;

    X.fill(0);
    int O = 1;
    std::array<BinaryInt, order+order_plus+1> orders_num;
    for( long t_MC = 0; t_MC<N_iter; ++t_MC){
      auto const P = compute_P(X);
      auto [times,
	    X_new] = invoke_heat_bath_mc<order+order_plus, order, is_translation_invariant>(P,
											  Phi,
											  X,
											  O);



      std::array<Real, order+order_plus+1> orders_temp; orders_temp.fill(0.);
      orders_num.fill(0);
      for(BinaryInt S=0; S<(1<<(order+order_plus)); ++S){
	auto const card_S = __builtin_popcount(S);
	orders_num[card_S]++;
	orders_temp[card_S] += times[S];
      }
      for(int j=0; j <= order+order_plus; ++j){
	orders[j] += orders_num[order]*orders_temp[j]/orders_num[j];
      }

      
      X = X_new;
    }
    std::cerr<<"order = "<<order<<", order_plus = "<<order_plus<<", L = "<<num_sites<<std::endl;
    auto [normalization, err_norm] = orders[0].MeanValueError();
    if constexpr(is_translation_invariant){
	normalization /= num_sites;
	err_norm /= num_sites;
      }
    std::array<Real, order+order_plus+1> variances, autocorr;
    for(int j=0; j <= order+order_plus; ++j){
      using ffd::core_math::Factorial;
      autocorr[j] = orders[j].AutocorrelationTime();
      auto [mean, error] = orders[j].MeanValueError();
      variances[j] = orders[j].Variance()/normalization/normalization;
      error += std::abs(mean*(err_norm/normalization));
      std::cerr<<j<<": "<<mean/normalization<<"+/-"<<error/normalization
	       <<", \u03C4 = "<< autocorr[j]
	       <<", \u03C3\u00B2 = "<<variances[j]<<", \u03C3\u00B2_eff = "<<variances[j]*(1+2*autocorr[j])<<std::endl;
    }

  
    {
      int constexpr j = 1;
      std::array<int, j> L_vec;
      L_vec.fill(num_sites);
      Real exact_f = 0.;
      Real exact_af = 0.;
      for( auto r: ffd::user_space::CartesianProductRange(L_vec) ){
	auto f_ex = f(r);
	exact_f += f_ex;
	exact_af += std::abs(f_ex);
      }
      std::cerr<<j<<": "<<exact_f<<" "<<exact_af/exact_f<<" "<<(exact_af*exact_af-exact_f*exact_f)/variances[j]/orders_num[j]<<" "<<orders_num[j]<<std::endl;
    }


        {
      int constexpr j = 2;
      std::array<int, j> L_vec;
      L_vec.fill(num_sites);
      Real exact_f = 0.;
      Real exact_af = 0.;
      for( auto r: ffd::user_space::CartesianProductRange(L_vec) ){
	auto f_ex = f(r);
	exact_f += f_ex;
	exact_af += std::abs(f_ex);
      }
      std::cerr<<j<<": "<<exact_f<<" "<<exact_af/exact_f<<" "<<(exact_af*exact_af-exact_f*exact_f)/variances[j]/orders_num[j]<<" "<<orders_num[j]<<std::endl;
    }



	        {
      int constexpr j = 3;
      std::array<int, j> L_vec;
      L_vec.fill(num_sites);
      Real exact_f = 0.;
      Real exact_af = 0.;
      for( auto r: ffd::user_space::CartesianProductRange(L_vec) ){
	auto f_ex = f(r);
	exact_f += f_ex;
	exact_af += std::abs(f_ex);
      }
      std::cerr<<j<<": "<<exact_f<<" "<<exact_af/exact_f<<" "<<(exact_af*exact_af-exact_f*exact_f)/variances[j]/orders_num[j]<<" "<<orders_num[j]<<std::endl;
    }


        {
      int constexpr j = 4;
      std::array<int, j> L_vec;
      L_vec.fill(num_sites);
      Real exact_f = 0.;
      Real exact_af = 0.;
      for( auto r: ffd::user_space::CartesianProductRange(L_vec) ){
	auto f_ex = f(r);
	exact_f += f_ex;
	exact_af += std::abs(f_ex);
      }
      std::cerr<<j<<": "<<exact_f<<" "<<exact_af/exact_f<<" "<<(exact_af*exact_af-exact_f*exact_f)/variances[j]/orders_num[j]<<" "<<orders_num[j]<<std::endl;
    }

	        {
		  int constexpr j = 5;
      std::array<int, j> L_vec;
      L_vec.fill(num_sites);
      Real exact_f = 0.;
      Real exact_af = 0.;
      for( auto r: ffd::user_space::CartesianProductRange(L_vec) ){
	auto f_ex = f(r);
	exact_f += f_ex;
	exact_af += std::abs(f_ex);
      }
      std::cerr<<j<<": "<<exact_f<<" "<<exact_af/exact_f<<" "<<(exact_af*exact_af-exact_f*exact_f)/variances[j]/orders_num[j]<<" "<<orders_num[j]<<std::endl;
    }

		
		{
		  int constexpr j = 6;
      std::array<int, j> L_vec;
      L_vec.fill(num_sites);
      Real exact_f = 0.;
      Real exact_af = 0.;
      for( auto r: ffd::user_space::CartesianProductRange(L_vec) ){
	auto f_ex = f(r);
	exact_f += f_ex;
	exact_af += std::abs(f_ex);
      }
      std::cerr<<j<<": "<<exact_f<<" "<<exact_af/exact_f<<" "<<(exact_af*exact_af-exact_f*exact_f)/variances[j]/orders_num[j]<<" "<<orders_num[j]<<std::endl;
    }


				{
		  int constexpr j = 7;
      std::array<int, j> L_vec;
      L_vec.fill(num_sites);
      Real exact_f = 0.;
      Real exact_af = 0.;
      for( auto r: ffd::user_space::CartesianProductRange(L_vec) ){
	auto f_ex = f(r);
	exact_f += f_ex;
	exact_af += std::abs(f_ex);
      }
      std::cerr<<j<<": "<<exact_f<<" "<<exact_af/exact_f<<" "<<(exact_af*exact_af-exact_f*exact_f)/variances[j]/orders_num[j]<<" "<<orders_num[j]<<std::endl;
    }

				
				{
				  int constexpr j = 8;
      std::array<int, j> L_vec;
      L_vec.fill(num_sites);
      Real exact_f = 0.;
      Real exact_af = 0.;
      for( auto r: ffd::user_space::CartesianProductRange(L_vec) ){
	auto f_ex = f(r);
	exact_f += f_ex;
	exact_af += std::abs(f_ex);
      }
      std::cerr<<j<<": "<<exact_f<<" "<<exact_af/exact_f<<" "<<(exact_af*exact_af-exact_f*exact_f)/variances[j]/orders_num[j]<<" "<<orders_num[j]<<std::endl;
    }


    // 								{
    // 				  int constexpr j = 9;
    //   std::array<int, j> L_vec;
    //   L_vec.fill(num_sites);
    //   Real exact_f = 0.;
    //   Real exact_af = 0.;
    //   for( auto r: ffd::user_space::CartesianProductRange(L_vec) ){
    // 	auto f_ex = f(r);
    // 	exact_f += f_ex;
    // 	exact_af += std::abs(f_ex);
    //   }
    //   std::cerr<<j<<": "<<exact_f<<" "<<exact_af/exact_f<<" "<<(exact_af*exact_af-exact_f*exact_f)/variances[j]/orders_num[j]<<" "<<orders_num[j]<<std::endl;
    // }


    // 				{
    // 				  int constexpr j = 10;
    //   std::array<int, j> L_vec;
    //   L_vec.fill(num_sites);
    //   Real exact_f = 0.;
    //   Real exact_af = 0.;
    //   for( auto r: ffd::user_space::CartesianProductRange(L_vec) ){
    // 	auto f_ex = f(r);
    // 	exact_f += f_ex;
    // 	exact_af += std::abs(f_ex);
    //   }
    //   std::cerr<<j<<": "<<exact_f<<" "<<(exact_af*exact_af-exact_f*exact_f)/variances[j]/orders_num[j]<<" "<<orders_num[j]<<std::endl;
    // }


    // 								{
    // 				  int constexpr j = 11;
    //   std::array<int, j> L_vec;
    //   L_vec.fill(num_sites);
    //   Real exact_f = 0.;
    //   Real exact_af = 0.;
    //   for( auto r: ffd::user_space::CartesianProductRange(L_vec) ){
    // 	auto f_ex = f(r);
    // 	exact_f += f_ex;
    // 	exact_af += std::abs(f_ex);
    //   }
    //   std::cerr<<j<<": "<<exact_f<<" "<<(exact_af*exact_af-exact_f*exact_f)/variances[j]/orders_num[j]<<" "<<orders_num[j]<<std::endl;
    // }

    // 				{
    // 				  int constexpr j = 12;
    //   std::array<int, j> L_vec;
    //   L_vec.fill(num_sites);
    //   Real exact_f = 0.;
    //   Real exact_af = 0.;
    //   for( auto r: ffd::user_space::CartesianProductRange(L_vec) ){
    // 	auto f_ex = f(r);
    // 	exact_f += f_ex;
    // 	exact_af += std::abs(f_ex);
    //   }
    //   std::cerr<<j<<": "<<exact_f<<" "<<(exact_af*exact_af-exact_f*exact_f)/variances[j]/orders_num[j]<<" "<<orders_num[j]<<std::endl;
    // }
    // 				{
    // 				  int constexpr j = 13;
    //   std::array<int, j> L_vec;
    //   L_vec.fill(num_sites);
    //   Real exact_f = 0.;
    //   Real exact_af = 0.;
    //   for( auto r: ffd::user_space::CartesianProductRange(L_vec) ){
    // 	auto f_ex = f(r);
    // 	exact_f += f_ex;
    // 	exact_af += std::abs(f_ex);
    //   }
    //   std::cerr<<j<<": "<<exact_f<<" "<<(exact_af*exact_af-exact_f*exact_f)/variances[j]/orders_num[j]<<" "<<orders_num[j]<<std::endl;
    // }

  }


}//namespace
