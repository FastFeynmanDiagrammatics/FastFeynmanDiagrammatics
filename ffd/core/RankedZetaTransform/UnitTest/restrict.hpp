namespace ffd::ranked_zeta::unit_test{

	void test_restrict(){
		using NilpotentPolynomial = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Full>;
		// using EvenNilPoly = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Even>;
		using namespace user_space;
		using std::size;

		const int order = 10;

		NilpotentPolynomial P1(order);
		for(int j=0; j < 1<<order; ++j){
			P1[j] = std::sin(0.7+4.4*j);
		}

		auto P1_r = P1.restrict_to_set(20);

		RankedZeta Z1 = P1;

		auto Z1_r = Z1.restrict_to_set(20);

		NilpotentPolynomial P1_r_new = Z1_r;

		for(BinaryInt j=0; j < P1_r.size(); ++j){
			// std::cerr << j << " " <<  P1_r[j] << " " << P1_r_new[j] << " " << std::abs(P1_r[j] - P1_r_new[j]) << std::endl;
			if (P1_r[j]!=0. || P1_r_new[j]!=0.) assert(std::abs(P1_r[j] - P1_r_new[j]) <= std::max(P1_r[j] * 1.e-10, 1.e-10) );
		}
	}

}//namespace
