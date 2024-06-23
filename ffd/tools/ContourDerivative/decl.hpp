namespace ffd::contour_derivative{

  namespace parameters{
    Real const RadiusPoleDefault = .1;
    unsigned long const NumberStartingPointsIntegrationDefault = 1ul<<10;
    Real const AbsolutePrecisionDefault = 1e-12;
  }
  

  
  template<typename complex_function_t>
  struct ContourDerivative{
    complex_function_t ComplexFunction;
    Real RadiusPole, AbsolutePrecision;
    unsigned long NumberStartingPointsIntegration =
      parameters::NumberStartingPointsIntegrationDefault;


    ContourDerivative(complex_function_t ComplexFunction_,
		      Real RadiusPole_ = parameters::RadiusPoleDefault,
		      Real AbsolutePrecision_ = parameters::AbsolutePrecisionDefault):
      ComplexFunction(ComplexFunction_), RadiusPole(RadiusPole_),
      AbsolutePrecision(AbsolutePrecision_) {}

    
    
    Complex operator()(int order,
		       Complex z0 = Complex(0., 0.)){
      using ffd::core_math::Pi;
      using ffd::user_space::expI;
      using ffd::vector_range::Range;
      using ffd::binary_accumulator::BinaryAccumulator;
      
      
      auto NumberPointsIntegration = NumberStartingPointsIntegration;
      Complex der_old = 0., der_now = 0.;
      do{
	der_old = der_now;
	der_now = 0;
	NumberPointsIntegration *= 2;
	BinaryAccumulator<Complex> Der_now_accumulator;
	for( auto j: Range(NumberPointsIntegration) ){
	  Complex zOnTheCircle = RadiusPole*expI(2*Pi*(j+0.5)/NumberPointsIntegration);
	  Der_now_accumulator += pow(zOnTheCircle, -order)*ComplexFunction(z0+zOnTheCircle);
	}
	der_now = Der_now_accumulator.MeanValue();
      }while( false
	      //std::abs(der_old-der_now) < AbsolutePrecision
	      );

      
      return der_now*ffd::core_math::Factorial(order);
    }


  };

  

  template<typename complex_function_t>
  Complex Derivative(complex_function_t F,
		     int order,
		     Complex z0 = Complex(0., 0.),
		     Real RadiusPole = parameters::RadiusPoleDefault,
		     Real AbsolutePrecision = parameters::AbsolutePrecisionDefault){
    ContourDerivative Der(F, RadiusPole, AbsolutePrecision);
    return Der(order, z0);
  }
  
	
}//namespace
