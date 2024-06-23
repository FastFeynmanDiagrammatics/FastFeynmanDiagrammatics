namespace ffd::ranked_zeta::unit_test{

void test_shift(){
	using NilpotentPolynomial = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Full>;
	// using EvenNilPoly = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Even>;
	using namespace user_space;
	using std::size;

	const int order = 10;

	// NilpotentPolynomial P(order);
	// P[20] = 1.;
	//
	// RankedZeta Z = P;

  //
  // for (BinaryInt set=0; set<Z.cardinality; ++set)
  // {
	//   std::cerr << set << "   ";
	//   for (BinaryInt k=0; k<=Z.rank_size; ++k)
	//   {
	// 	  std::cerr << Z[set+Z.cardinality*k] << " ";
	//   }
	//   std::cerr << std::endl;
  // }
  // exit(1);


	NilpotentPolynomial P1(order), P2(order);
	int set1 = 20;
	int set2 = 139;
	for(int j=0; j < 1<<order; ++j){
		if ((set1 & j)==j){
			P1[j] = std::sin(0.7+4.4*j);
		}
		if ((set2 & j)==j){
			P2[j] = std::sin(1.5+2.3*j);
			// std::cerr << j << " " << P2[j] << std::endl;
		}
	}

	auto P1_sh = P1.shift(9,10);
	// auto P2_sh = P2.shift(2816,12);

	RankedZeta Z1 = P1;
	RankedZeta Z2 = P2;

	auto Z1_sh = Z1.shift(9,10);
	// auto Z2_sh = Z2.shift(2816,12);

	NilpotentPolynomial P1_sh_new = Z1_sh;
	// NilpotentPolynomial P2_sh_new = Z2_sh;

	for(int j=0; j < (1<<10); ++j){
		// std::cerr << j << " " <<  P1_sh[j] << " " << P1_sh_new[j] << " " << std::abs(P1_sh[j] - P1_sh_new[j]) << std::endl;
		if (P1_sh[j]!=0. || P1_sh_new[j]!=0.) assert(std::abs(P1_sh[j] - P1_sh_new[j]) <= std::max(P1_sh[j] * 1.e-10, 1.e-10) );
	}
	// for(int j=0; j < (1<<12); ++j){
	// 	// std::cerr << j << " " <<  P2_sh[j] << " " << P2_sh_new[j] << " " << std::abs(P2_sh[j] - P2_sh_new[j]) << std::endl;
	// 	if (P2_sh[j]!=0. || P2_sh_new[j]!=0.) assert(std::abs(P2_sh[j] - P2_sh_new[j]) <= std::max(P2_sh[j] * 1.e-10, 1.e-10) );
	// }
}

}//namespace
