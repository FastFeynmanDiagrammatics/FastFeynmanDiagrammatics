namespace ffd::determinant::unit_test{


  void matrix3(){

    std::vector<Real> Matrix3{0.958851077208406, 1.788435374475382, 2.2333204859216336, 2.532414720122871, 
 2.727116370834386, 2.828323306541442, 2.8404724780477943, 2.767600175374829, 
			      2.615307622186636};


    Real det_ex = 0.2089277028761125;


    auto det_nm = Determinant(Matrix3);


    assert((
	    std::abs(det_ex-det_nm) <
	    10*std::numeric_limits<Real>::epsilon()
	    ));
    
  }
  

}//namespace
