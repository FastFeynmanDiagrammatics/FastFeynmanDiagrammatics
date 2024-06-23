

namespace ffd::user_space::imaginary_time_lattice_propagator{

  std::vector<std::vector<Complex>>
  HalveMatrix(std::vector<std::vector<Complex>> const& H){
    using ffd::vector_range::Range;
    
    int const size_half_H = std::size(H)/2;
    std::vector<std::vector<Complex>>
      half_H(size_half_H,
	     std::vector<Complex>(size_half_H, 0.));

    for( int j: Range(size_half_H) ){
      for( int m: Range(size_half_H) ){
	half_H[j][m] = H[j][m];
      }
    }

    return half_H;
  }


}//namespace
