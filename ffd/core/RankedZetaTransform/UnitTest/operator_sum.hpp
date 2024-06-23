#pragma once

namespace ffd::ranked_zeta::unit_test{

  void operator_sum(){

	  using NilpotentPolynomial = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Full>;
	  using EvenNilPoly = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Even>;
	  using namespace user_space;

	  const int order = 10;
	  NilpotentPolynomial P1(order), P2(order);
	  for(int j=0; j < 1<<order; ++j){
		  P1[j] = std::sin(0.7+4.4*j);
		  P2[j] = std::sin(1.5+2.3*j);
		}

	auto P_sum = P1+P2;
	auto P_dif = P1-P2;

	RankedZeta Z1 = P1;
	RankedZeta Z2 = P2;

	auto Z_sum = Z1+Z2;
	auto Z_dif = Z1-Z2;

	NilpotentPolynomial P_sum_new = Z_sum;
	NilpotentPolynomial P_dif_new = Z_dif;

	for(int j=0; j < 1<<order; ++j){
	   assert(std::abs(P_sum[j] - P_sum_new[j]) < 100*std::numeric_limits<long double>::epsilon() );
	   assert(std::abs(P_dif[j] - P_dif_new[j]) < 100*std::numeric_limits<long double>::epsilon() );
	}

	EvenNilPoly P1_e(order), P2_e(order);
	for(int j=0; j < 1<<order; ++j){
		if (!__builtin_parity(j)){
			P1_e[j] = std::sin(0.7+4.4*j);
			P2_e[j] = std::sin(1.5+2.3*j);
		}
	  }

  auto P_sum_e = P1_e+P2_e;
  auto P_dif_e = P1_e-P2_e;

  EvenZeta Z1_e = P1_e;
  EvenZeta Z2_e = P2_e;

  auto Z_sum_e = Z1_e+Z2_e;
  auto Z_dif_e = Z1_e-Z2_e;

  EvenNilPoly P_sum_new_e = Z_sum_e;
  EvenNilPoly P_dif_new_e = Z_dif_e;

  for(int j=0; j < 1<<order; ++j){
	 assert(std::abs(P_sum_e[j] - P_sum_new_e[j]) < 100*std::numeric_limits<long double>::epsilon() );
	 assert(std::abs(P_dif_e[j] - P_dif_new_e[j]) < 100*std::numeric_limits<long double>::epsilon() );
  }

  }

}//namespace
