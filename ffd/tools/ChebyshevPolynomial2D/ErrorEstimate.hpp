

namespace ffd::chebyshev_polynomial_2d{

  template<typename Field>
  Real ChebyshevPolynomial2D<Field>::ErrorEstimate() const{
    Real MaximalAbsoluteError = 0;

    using namespace std;
    
    if(Order == 1){
      return MaximalAbsoluteError = abs(Coef[0]);
    }

    for(int position=0; position < Order; ++position){
      for(bool x_direction: {false, true}){
	int kx   = x_direction ? position : Order-1;
	int ky   = x_direction ? Order-1  : position;
	int kx_m = x_direction ? position : Order-2;
	int ky_m = x_direction ? Order-2  : position;
	MaximalAbsoluteError = max({MaximalAbsoluteError,
				    abs(Coef[kx+Order*ky]),
				    abs(Coef[kx_m+Order*ky]),
				    abs(Coef[kx_m+Order*ky_m]),
				    abs(Coef[kx+Order*ky_m])});
      }
    }
    return MaximalAbsoluteError;
  }


}//namespace
