#ifndef FFD_COMPLEX_ZEROS_HEADER_NOT_BEEN
#define FFD_COMPLEX_ZEROS_HEADER_NOT_BEEN


namespace ffd::ComplexZeros{
  
  struct ComplexZeros{

    
    std::function<Complex(Complex)> Function;

    
    ffd::AnalyticalDerivative::AnalyticalDerivative Derivative;


    
    Real PrecisionOfZero, EpsilonStep = 0.05;

    
    ComplexZeros(std::function<Complex(Complex)> func, Real radius_der = 0.001, Real Epsilon_Step = 0.01, Real precision_zero = 1e-13): Function(func), Derivative(func, radius_der), PrecisionOfZero(precision_zero), EpsilonStep(Epsilon_Step) {}

    
    Complex FindNearestZero(){
      Complex z0{1e-30, 1e-30};
      while(std::abs(Function(z0)) > PrecisionOfZero){
	Complex Der1 = Derivative.CoefficientOfzTo1(1, z0);
	Complex Der2 = Derivative.CoefficientOfzTo1(2, z0);
	Complex ValueFunc = Function(z0);
	Complex Step = (ValueFunc/Der1)*EpsilonStep*(-1.);
	if(std::abs(Step) > 10){
	  Step = std::sqrt(ValueFunc/Der2)*EpsilonStep*(-1.);
	}
	z0 += Step;
      }
      return z0;
    }
    
  };
  
}//namespace ffd::ComplexZeros

#endif
