

#include"returns_true_operator_parentheses.hpp"

namespace ffd::chebyshev_polynomial_2d::unit_test{

  void operator_parentheses(){
    using namespace std;

    Real xm = 0, xM = 1, ym = 0, yM = 1;
    Real precision = 1e-10;


    

    auto f = [](Real x, Real y){return (x+y*y*x*x)*tanh(x*y)*sin(x*y*y);};
    
    assert(returns_true_operator_parentheses(f,
					     {{{xm, xM}, {ym, yM}}},
					     precision));


    

    
    auto Sin_not_so_easy = [](Real x, Real y)->Real{return tanh(2*sin(2*x)-cos(y));};
    
    assert(returns_true_operator_parentheses(Sin_not_so_easy,
					     {{{-1, 1}, {-2, 2}}},
					     precision=1e-10));

    

    //compile with -Ofast to get this one
    // auto Sin_diabulus = [](Real x, Real y)->Real{return sin(13*sin(11*x)-12*cos(13*y))*tanh(2*x*y+x*x-y+1);};
    // assert(returns_true_operator_parentheses(Sin_diabulus,
    // 					     {{{0, 1}, {-1, 0}}},
    // 					     precision=1e-12));



  }

}//namespace
