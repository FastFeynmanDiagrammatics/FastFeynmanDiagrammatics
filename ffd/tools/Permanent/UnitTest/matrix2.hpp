namespace ffd::permanent::unit_test{

  void matrix2(){

    std::array<Real, 4> Matrix2{0.958851077208406,
				1.788435374475382,
				2.2333204859216336,
				2.532414720122871};

    Real per_ex = 6.422357941891038;

    
    auto per_nm = Permanent(Matrix2);


    //std::cerr<<per_nm<<" "<<per_ex-per_nm<<std::endl;

    assert((
	    std::abs(per_ex-per_nm) < 1e-10
	    ));

  }


}//namespace
