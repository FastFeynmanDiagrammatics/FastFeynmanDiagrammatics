

#include"easy_to_fit.hpp"
#include"not_so_easy_to_fit.hpp"
#include"harder_than_not_so_easy.hpp"
#include"compute_vector_of_values.hpp"
#include"compute_max_error.hpp"

namespace ffd::chebyshev_polynomial::unit_test{
  
  void operator_parentheses(){
    using std::size;
    const int N_Cheby = 40;
    easy_to_fit func(1.);
    ChebyshevPolynomial<Real> P(func, N_Cheby);

    auto max_diff = compute_max_error(P, func);
    assert( max_diff < N_Cheby*std::numeric_limits<Real>::epsilon() );
    assert( std::abs(P.ErrorEstimate()) < N_Cheby*std::numeric_limits<Real>::epsilon() );


    const Real epsilon = 1e-8;
    ChebyshevPolynomial<Real> P_auto(func, {-1, 1}, epsilon);

    auto max_diff_auto = compute_max_error(P_auto, func);
    assert( max_diff_auto <  epsilon);

    
    const int N_Cheby_not_so_easy = 350;
    not_so_easy_to_fit SinOfSin(10., 10.);
    ChebyshevPolynomial<Real> P2(SinOfSin, N_Cheby_not_so_easy);

    auto max_diff2 = compute_max_error(P2, SinOfSin);
    assert( max_diff2 < N_Cheby_not_so_easy*std::numeric_limits<Real>::epsilon() );
    assert( std::abs(P2.ErrorEstimate()) < N_Cheby_not_so_easy*std::numeric_limits<Real>::epsilon() );

    
    const Real epsilon2 = 1e-8;
    ChebyshevPolynomial<Real> P2_auto(SinOfSin, {-1, 1}, epsilon2);

    auto max_diff2_auto = compute_max_error(P2_auto, SinOfSin);
    assert( max_diff2_auto <  100*epsilon2);

    
    const int N_Cheby_harder = 900;
    harder_than_not_so_easy SinOfSinOfSinTimesGaussianExp{};  //impossibly-looking function
    ChebyshevPolynomial<Real> P3(SinOfSinOfSinTimesGaussianExp, N_Cheby_harder);

    auto max_diff3 = compute_max_error(P3, SinOfSinOfSinTimesGaussianExp);
    assert( max_diff3 < N_Cheby_harder*std::numeric_limits<Real>::epsilon() );
    assert( std::abs(P3.ErrorEstimate()) < N_Cheby_harder*std::numeric_limits<Real>::epsilon() );


    const Real epsilon3 = 1e-8;
    ChebyshevPolynomial<Real> P3_auto(SinOfSinOfSinTimesGaussianExp, {-1, 1}, epsilon3);

    auto max_diff3_auto = compute_max_error(P3_auto, SinOfSinOfSinTimesGaussianExp);
    assert( max_diff3_auto <  100*epsilon3);


  }

}//namespace
