namespace ffd::ranked_zeta::unit_test{

	void scalar_operations(){
		using NilpotentPolynomial = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Full>;
		using namespace user_space;

	  const int order = 10;
	  NilpotentPolynomial P(order);
	  NilpotentPolynomial Q_;
	  for(int j=0; j < 1<<order; ++j){
		P[j] = j;
	  }

	  RankedZeta Z = P;

	  auto Q = P;
	  auto R = Z;
	  long double x = 1.23456789238123876213l;

	  Q += x;
	  R += x;
	  Q_ = R;
	  for(int j=0; j < 1<<order; ++j){
		assert(std::abs(Q[j] - Q_[j]) < 10*std::numeric_limits<Real>::epsilon());
	  }

	  Q = P - x;
  	  R = Z - x;
	  Q_ = R;
	  for(int j=0; j < 1<<order; ++j){
		assert(std::abs(Q[j] - Q_[j]) < 10*std::numeric_limits<Real>::epsilon());
	  }

	  Q = P;
	  Q *= x;
	  R = Z;
	  R *= x;
	  Q_ = R;
	  for(int j=0; j < 1<<order; ++j){
		assert(std::abs(Q[j] - Q_[j]) < 10*(std::numeric_limits<Real>::min()+std::abs(Q[j])*std::numeric_limits<Real>::epsilon()));
	  }

	  Q = P/x;
	  R = Z/x;
	  Q_ = R;
	  for(int j=0; j < 1<<order; ++j){
		assert(std::abs(Q[j] - Q_[j]) < 10*(std::numeric_limits<Real>::min()+std::abs(Q[j])*std::numeric_limits<Real>::epsilon()));
	  }

	  Q = -P;
	  R = -Z;
	  Q_ = R;
	  for(int j=0; j < 1<<order; ++j){
		assert(std::abs(Q[j] - Q_[j]) < 10*std::numeric_limits<Real>::epsilon());
	  }

	}

}//namespace
