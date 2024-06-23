namespace ffd::analytical_derivative{
  //This class can be used to compute coefficients of
  //Taylor expansion of a function for which we have
  //an explicit analytic expression (to be stored in `value`)
  struct AnalyticalDerivative{
    std::function<Complex(Complex)> value;

    AnalyticalDerivative(){}

    AnalyticalDerivative(std::function<Complex(Complex)> func,
			 Real RadiusPoles1 = 0.1): value(func), RadiusPoles(RadiusPoles1) {}
    
    Real RadiusPoles = 0.1;

    long NumberPointsIntegration = (1<<14);


    Complex CoefficientOfzTo1(int order,
			      Complex z0 = Complex(0., 0.),
			      Real Radius = -1.){
      if(Radius < 0){
	Radius = RadiusPoles;
      }
      Complex ret = 0;
      for(long j=0; j<NumberPointsIntegration; ++j){
	Complex zOnTheCircle = Radius*exp(Complex(0, 2*ffd::core_math::Pi*(j+0.5)/NumberPointsIntegration));
	ret += pow(zOnTheCircle, -J)*value(z0+zOnTheCircle);
      }
      return ret /= NumberPointsIntegration;
    }

  };

}//namespace ffd::analytical_derivative
