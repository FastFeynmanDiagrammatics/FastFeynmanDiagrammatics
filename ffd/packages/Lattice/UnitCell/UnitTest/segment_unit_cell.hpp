

namespace ffd::lattice::unit_test{

  void segment_unit_cell(){

    SegmentUnitCell<3> S(1.5);

    std::array<std::array<Real, 3>, 2> offsets;
    offsets[0] = {0., 0., 0.};
    offsets[1] = {1.5, 0., 0.};
    for(auto j: {0, 1}){
      for(auto k: {0, 1, 2}){
	assert( std::abs(S.AtomicOffsets[j][k] - offsets[j][k]) <
		std::numeric_limits<Real>::epsilon() );
      }
    }
    
  }

}//namespace
