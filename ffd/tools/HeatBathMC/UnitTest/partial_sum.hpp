namespace ffd::heat_bath_mc::unit_test{

  void partial_sum_test(){
	  
    std::array<Real,(1<<4)> P;
    for(BinaryInt S=0; S<(1<<4); ++S){
      P[S]=1.;
    }

    auto R = ffd::heat_bath_mc::partial_sum<4,4>(P,1.);

    assert(R[0] == 1);
    for(BinaryInt S=1; S<(1<<4); ++S){
      //std::cout << S << " " << P[S] << " " << R[S] << " " << __builtin_popcount(S) << std::endl;
      assert(R[S] - __builtin_popcount(S) < 1.e-10);
    }
  }

}//namespace
