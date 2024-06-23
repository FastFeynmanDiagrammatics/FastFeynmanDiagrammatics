

namespace ffd::chebyshev_polynomial_2d::unit_test{

  bool returns_true_operator_parentheses(std::function<Real(Real, Real)> f,
					 std::array<std::array<Real, 2>, 2> interval,
					 Real precision){
    bool IsOk = true;
    

    using namespace std;


    ChebyshevPolynomial2D<Real> P(f,
				  interval,
				  precision);

    
    const int N_samples = 10;
    for(int x=0; x<N_samples; ++x){
      for(int y=0; y<N_samples; ++y){
	std::array<Real, 2> X;
	for(int j: {0, 1}){
	  X[j] = interval[j][0]
	    + ( (j==0)*x + (j==1)*y +.11231231)*(interval[j][1]-interval[j][0])/N_samples;
	}
	IsOk = IsOk && abs(f(X[0], X[1]) - P(X[0], X[1])) < precision;
	//	std::cerr<<"("<<X[0]<<", "<<X[1]<<") "<<f(X[0], X[1])<<" "<<(f(X[0], X[1]) - P(X[0], X[1]))<<std::endl;
      }
    }
    //std::cerr<<"order = "<<P.Order<<std::endl;
    
    return IsOk;

  }


}//namespace

