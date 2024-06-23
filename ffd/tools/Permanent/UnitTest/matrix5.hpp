namespace ffd::permanent::unit_test{


  void matrix5(){

    std::array<Real, 25> Matrix5{0.958851077208406, 1.788435374475382, 2.2333204859216336, 2.532414720122871, 
 2.727116370834386, 2.828323306541442, 2.8404724780477943, 2.767600175374829, 
 2.615307622186636, 2.3914104243534404, 2.106035197298822, 1.771426427134326, 
 1.4015755815048871, 1.0117327534380096, 0.617841945046836, 
 0.23593354462076033, -0.11849581122875152, -0.4310878739235038, 
 -0.689185801076189, -0.8823318734989105, -1.0026792829492417, 
 -1.045301060582747, -1.0083834861182301, -0.8932960313221309, 
				 -0.7045348844478017};




    Real per_ex = 351.5849076106582;


    auto per_nm = Permanent(Matrix5);

    //std::cerr<<per_nm<<" "<<per_ex-per_nm<<std::endl;

    assert((
	    std::abs(per_ex-per_nm) < 1e-10
	    ));
    
  }
  

}//namespace
