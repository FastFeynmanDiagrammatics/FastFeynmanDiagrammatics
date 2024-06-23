

namespace ffd::user_space::imaginary_time_lattice_propagator{

  template<int d>

  Real

  ScalarProductArrays(std::array<Real, d> x,
		     std::array<Real, d> y){
    Real scalar_product = 0;


    for( int j: ffd::vector_range::Range(d) ){
      scalar_product += x[j]*y[j];
    }


    return scalar_product;
    
  }

}//namespace
