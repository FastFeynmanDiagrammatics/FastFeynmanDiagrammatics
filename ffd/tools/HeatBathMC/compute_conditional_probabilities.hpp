namespace ffd::heat_bath_mc{

  template<int max_order,
	   bool is_translation_invariant>

  std::array<Real, (1<<max_order)>

  compute_conditional_probabilities(std::array<Real,
				    (max_order+!is_translation_invariant)*
				    (max_order+!is_translation_invariant)> const& seeds){
    using namespace ffd::set_theory;
    int constexpr eff_order = max_order+!is_translation_invariant;
    // BinaryInt constexpr two_eff_order = (1<<eff_order);
    BinaryInt constexpr two_max_order = (1<<max_order);
    std::array<Real, two_max_order> P_cond;
    P_cond.fill(0.);
    P_cond[0] = 1.;

    // timer_MC_6.ini();

    // for(BinaryInt S=1; S<two_max_order-is_translation_invariant; ++S){
    //   BinaryInt const card_S = __builtin_popcount(S);
    //   BinaryInt const card_cS = eff_order - card_S;
    //   std::vector<BinaryInt> const digits_S = VectorOfBinaryDigitsOf(S);
    //   std::vector<BinaryInt> const digits_cS =
    //VectorOfBinaryDigitsOf(two_eff_order-1-S);
    //   Real sum_cumul = 0.;
    //   for(BinaryInt s=0; s<card_S; ++s){
    // 	BinaryInt const s_i = digits_S[s];
    // 	BinaryInt const s_eff_index = s_i*eff_order;
    // 	Real pseudo_potential = 0.;
    // 	for(BinaryInt c=0; c<card_cS; ++c){
    // 	  pseudo_potential += seeds[s_eff_index+digits_cS[c]];
    // 	}
    // 	sum_cumul += P_cond[S-(1<<s_i)]*pseudo_potential;
    //   }
    //   P_cond[S] = sum_cumul/(card_S*card_cS);
    // }

    // for(BinaryInt S=1; S<two_max_order-is_translation_invariant; ++S){
    //   Real sum_cumul = 0.;
    //   for(BinaryInt j=0; j<max_order; ++j){
    // 	BinaryInt const two_j = 1<<j;
    // 	if(S & two_j){
    // 	  BinaryInt const j_eff_index = j*eff_order;
    // 	  Real pseudo_potential = 0.;
    // 	  for(BinaryInt k=0; k<eff_order; ++k){
    // 	    if((S & (1<<k))==0){
    // 	      pseudo_potential += seeds[j_eff_index+k];
    // 	    }
    // 	  }
    // 	  sum_cumul += P_cond[S-two_j]*pseudo_potential;
    // 	}
    //   }
    //   BinaryInt const card_S = __builtin_popcount(S);
    //   BinaryInt const card_cS = eff_order - card_S;
    //   P_cond[S] = sum_cumul/(card_S*card_cS);
    // }


    std::vector<double> sum_factor(eff_order+1,0.);
    for (BinaryInt card=0; card<eff_order; ++card){
      sum_factor[card] = 1.0/(double)(card*(eff_order-card));
    }

    std::vector<double> high_factor(eff_order+1,0.);
    for (BinaryInt hs=0; hs<=eff_order; ++hs){
      BinaryInt const high_S_max = hs * eff_order;
      for(BinaryInt js=0; js<hs; ++js){
        high_factor[hs] += seeds[high_S_max+js];
      }
      for(BinaryInt js=hs+1; js<eff_order; ++js){
        high_factor[hs] += seeds[high_S_max+js];
      }
    }

    std::vector<double> pseudo_pot(two_max_order*eff_order,0.);

    for(BinaryInt S=1; S<two_max_order-is_translation_invariant; ++S){
      BinaryInt const high_S       = 31 - __builtin_clz(S);
      BinaryInt const two_high_S   = 1 << high_S;
      BinaryInt const S_ord        = S * eff_order;
      BinaryInt const S_hs         = S - two_high_S;
      BinaryInt const high_S_max   = high_S * eff_order;
      BinaryInt const S_ord_high_S = S_ord + high_S;
      BinaryInt const high_S_ord   = S_hs * eff_order;

      Real sum_cumul = 0.;

      for(BinaryInt S_ = S; S_!=two_high_S; ){
        BinaryInt const j = __builtin_ctz(S_);
        BinaryInt const two_j = 1<<j;
        BinaryInt const S_ord_j = S_ord+j;
        pseudo_pot[S_ord_j] = pseudo_pot[high_S_ord+j] - seeds[j*eff_order+high_S];
        sum_cumul += P_cond[S-two_j]*pseudo_pot[S_ord_j];
        S_ -= two_j;
      }

      if (S_hs==0){
        pseudo_pot[S_ord_high_S] = high_factor[high_S];
      }
      else{
        BinaryInt const sub_S = 31 - __builtin_clz(S_hs);
        pseudo_pot[S_ord_high_S] = pseudo_pot[(S-(1<<sub_S))*eff_order+high_S] - seeds[high_S_max+sub_S];
      }
      sum_cumul += P_cond[S_hs]*pseudo_pot[S_ord_high_S];

      P_cond[S] = sum_cumul*sum_factor[__builtin_popcount(S)];
    }

    // timer_MC_6.fin();

    if constexpr(is_translation_invariant){
	Real sum_cumul = 0.;
	for(int j=0; j<max_order; ++j){
	  sum_cumul += P_cond[two_max_order-1-(1<<j)];
	}
	P_cond[two_max_order-1] = sum_cumul/max_order;
      }


      return P_cond;
    }


  }//namespace
