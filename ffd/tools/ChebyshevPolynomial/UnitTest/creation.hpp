

namespace ffd::chebyshev_polynomial::unit_test{

  Real function_test(Real x){
    return x;
  }

  struct functor_test{
    Real a;
    functor_test(Real a_): a(a_) {}
    Real operator()(Real x){
      return a*x*x*x;
    }
  };

  
  void creation(){

    ChebyshevPolynomial P0{};

    ChebyshevPolynomial<Real> P1;

    static_assert(std::is_same<decltype(P0), decltype(P1)>::value);

    ChebyshevPolynomial<Complex> P2;

    std::vector<Real> vector_of_values = {1., 2. ,3. , 4.};
    
    ChebyshevPolynomial P3(vector_of_values);
    ChebyshevPolynomial P4(vector_of_values, 0, 1);

    ChebyshevPolynomial<Real> P5(function_test, 10);
    ChebyshevPolynomial<Real> P6(function_test, 10, -2, 3);

    auto lambda = [](Real x)->Real{return x*x-1;};
    ChebyshevPolynomial<Real> P7(lambda, 20);

    ChebyshevPolynomial<Real> P8(functor_test(-1.), 20);
    
  }

}//namespace
