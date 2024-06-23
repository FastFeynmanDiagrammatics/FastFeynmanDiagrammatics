namespace ffd::determinant::unit_test{

  void matrix2(){

    std::array<Real, 4> Matrix2{0.958851077208406,
				1.788435374475382,
				2.2333204859216336,
				2.532414720122871};

    Real det_ex = -1.5659407772345602;

    
    auto det_nm = Determinant(Matrix2);


    // std::cerr<<det_nm<<" "<<det_ex-det_nm<<std::endl;

    assert((
	    std::abs(det_ex-det_nm) < 10*std::numeric_limits<Real>::epsilon()
	    ));

  }


}//namespace
