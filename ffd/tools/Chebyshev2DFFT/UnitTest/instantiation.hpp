

namespace ffd::chebyshev_2d_fft::unit_test{

  Real func_test_simple(Real x, Real y){
    return exp(-2*x)*sin(10*y)*sin(x*x*y);
  }
  
  
  void
  instantiation(){
    Real precision = 1e-12;
    
    Chebyshev2DFFT P;
    Chebyshev2DFFT Q = 1.;

    std::array<std::array<Real, 2>, 2> bounds;
    bounds.fill( std::array{0., 1.} );

    Chebyshev2DFFT R(func_test_simple,
		     bounds,
		     precision);

    auto [x_size, y_size] = sizes(R);
    // std::cerr<<x_size<<" "<<y_size<<std::endl;
    // std::cerr<<std::setprecision(16)<<R(0.2, .4)<<" "<<func_test_simple(.2, .4)<<std::endl;
    
  }


}//namespace
