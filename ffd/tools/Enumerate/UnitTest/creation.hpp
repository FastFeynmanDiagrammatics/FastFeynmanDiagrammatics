namespace ffd::enumerate::unit_test{

  void creation(){

    Real z0 = 3.14, z1 = 6.28;
    
    std::array<Real, 2> x{z0, z1};
    
    auto y = Enumerate(x);

    int iter = 0;
    for( auto [j, z]: Enumerate(x) ){
      // std::cerr<<j<<" "<<z<<std::endl;

      assert((
	      j == iter
	      ));

      Real z_exact = z0 * (iter == 0) + z1 * (iter == 1);
	
      assert((
	      std::abs(z_exact - z) < std::numeric_limits<Real>::epsilon()
	      ));

      ++iter;
    }

  }

}//namespace
