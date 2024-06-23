namespace ffd::permanent::unit_test{


  void matrix3(){

    std::vector<Real> Matrix3{0.958851077208406, 1.788435374475382, 2.2333204859216336, 2.532414720122871, 
 2.727116370834386, 2.828323306541442, 2.8404724780477943, 2.767600175374829, 
			      2.615307622186636};


    Real per_ex = 73.50977075314128;


    auto per_nm = Permanent(Matrix3);
    
    //std::cerr<<per_nm<<" "<<per_ex-per_nm<<std::endl;
 
    assert((
	    std::abs(per_ex-per_nm) < 1e-10
	    ));
    
  }
  

}//namespace
