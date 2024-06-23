

namespace ffd::lattice::unit_test{

  void equilateral_triangle_unit_cell(){

    EquilateralTriangleUnitCell<4> T(1.6);

    std::array<std::array<Real, 4>, 3> offsets;
    offsets[0] = {0., 0., 0., 0.};
    offsets[1] = {1.6, 0., 0., 0.};
    offsets[2] = {.5*1.6, .5*std::sqrt(3.)*1.6, 0., 0.};
    for(int j=0; j < 3; ++j){
      for(int k=0; k < 4; ++k){
	assert( std::abs(offsets[j][k] - T.AtomicOffsets[j][k]) <
		std::numeric_limits<Real>::epsilon() );
      }
    }
    
  }
  
}//namespace
