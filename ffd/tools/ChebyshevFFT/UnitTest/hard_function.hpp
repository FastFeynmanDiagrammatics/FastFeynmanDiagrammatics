

namespace ffd::chebyshev_fft::unit_test{

  struct hard_function{
    std::array<Real, 5> omega{{11., 13., 7., 17., 19.}};
    std::array<Real, 5> shifts{{97., 29., 0., 37., 0.}};


    Real operator()(Real x) const{
      using std::sin, std::exp;
      return sin(shifts[0] +
		 omega[0]*sin(shifts[1] +
			      omega[1]*sin(shifts[2]+omega[2]*x)))*
	exp(exp(sin(shifts[3]+omega[3]*sin(omega[4]*x+shifts[4]))));
    }

  };


}//namespace
