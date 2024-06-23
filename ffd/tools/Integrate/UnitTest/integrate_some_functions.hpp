namespace ffd::integrate::unit_test{

  void integrate_some_functions(){
    using ffd::core_math::Pi;
    using ffd::core_math::Factorial;
    using std::pow, std::abs, std::exp, std::sin, std::cos, std::acos, std::log;

    
    Real epsilon = 1e-12;
    Real b = 2.9213;
    auto I = Integrate<Method::ClenshawCurtis>([](Real x){return exp(-x);}, {0, b}, epsilon);
    assert( abs( I - (1-exp(-b)) ) < epsilon );

    
    Real a = 3.243241;
    b = 2.9213;
    I = Integrate<Method::ClenshawCurtis>([a, b](Real x){return 1./(a+b*cos(x));}, {0, Pi/2}, epsilon);
    assert( abs( I - acos(b/a)/sqrt((a+b)*(a-b)) ) < epsilon );


    a = .32738127;
    I = Integrate<Method::ClenshawCurtis>([a](Real x){return x*sin(x)/(1-2*a*cos(x)+a*a);},
					  {0, Pi}, epsilon);
    assert( abs( I - Pi*log(1+a)/a ) < epsilon );


    a = .62738127;
    int m = 13;
    I = Integrate<Method::ClenshawCurtis>([a, m](Real x){return cos(m*x)/(1-2*a*cos(x)+a*a);},
					  {0, Pi}, epsilon);
    assert( abs( I - Pi*pow(a, m)/(1-a*a) ) < epsilon );

    
    m = 13;
    int n = 25;
    I = Integrate<Method::ClenshawCurtis>([m, n](Real x){return pow(x, m)*pow(log(x), n);},
					  {0, 1}, epsilon);
    assert( abs( I + pow(-1./(m+1), n+1)*Factorial(n)) < epsilon );

    
    m = 11;
    n = 9;
    I = Integrate<Method::ClenshawCurtis>([m, n](Real x){return pow(x, m)*pow(1-x, n);},
					  {0, 1}, epsilon);
    assert( abs( I -Factorial(n)*1./Factorial(m+n+1)*Factorial(m) ) < epsilon );

  }

}//namespace
