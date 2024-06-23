namespace ffd::heat_bath_mc::unit_test{

  void conformal_map_test(){
	  
    std::array<Real,(1<<4)> P;
    for(BinaryInt S=0; S<(1<<4); ++S){
      P[S]=1.;
    }

    std::vector<double> mat_conform(5*5,0.);
    for (BinaryInt y=0; y<5; ++y){
      for (BinaryInt x=0; x<=y; ++x){
        mat_conform[x+5*y]=1.;
      }
    }

    // ffd::math_tools::print_matrix(mat_conform);

    auto B = ffd::heat_bath_mc::conformal_map<4,4>(P, mat_conform);
    
    assert(B[0] == 1);
    for(BinaryInt S=1; S<(1<<4); ++S){
      // std::cout << S << " " << P[S] << " " << B[S] << " " << __builtin_popcount(S) << std::endl;
      assert(B[S] - __builtin_popcount(S) < 1.e-10);
    }
  }

}//namespace
