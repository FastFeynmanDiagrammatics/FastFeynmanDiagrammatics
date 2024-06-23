namespace ffd::gauss_pfaffian::unit_test{

  void four_by_four(){

    
    Real a = 1.21312,
      b = 2.21312,
      c = -.912312,
      d = 1.231212,
      e = -2.28822,
      f = 1.847327;
    
    
    std::array<Real, 6> mat{a, b, c, d, e, f};
    Real pfaffian_exact = a*f-b*e+c*d;


    assert((
	    std::abs(
		     pfaffian_exact -
		     Pfaffian(mat)
		     ) <
	    10*std::numeric_limits<Real>::epsilon()
	    ));
    
    
    // std::cerr<<Pfaffian(mat)<<std::endl;

    

    // std::cerr<<pfaffian_exact<<std::endl;

  }

}//namespace
