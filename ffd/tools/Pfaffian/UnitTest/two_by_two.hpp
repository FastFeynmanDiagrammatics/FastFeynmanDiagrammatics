namespace ffd::gauss_pfaffian::unit_test{

  void two_by_two(){

    
    std::vector<Real> matx;

    
    assert((
	    std::abs(Pfaffian(matx) - 1.) <
	    std::numeric_limits<Real>::epsilon()*10
	    ));

    
    std::array<Real, 1> mat{2.12312};

    
    // std::cerr<<Pfaffian(mat)<<std::endl;

    
    assert((
	    std::abs( mat[0] -
		      Pfaffian(mat)) <
	    std::numeric_limits<Real>::epsilon()*10
	    ));
    
  }


}//namespace
