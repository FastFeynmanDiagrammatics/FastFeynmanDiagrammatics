namespace ffd::periodic_coordinate::unit_test{

  void periodic_coordinates_real_space(){
    using ffd::merge, ffd::get, ffd::size, ffd::sign;
    
    PeriodicCoordinates<2, Real, int, int, int> X;
    get<1>(X).RealSpaceVectors = {{1., 0.}};
    get<2>(X).RealSpaceVectors = {{0., 1.}};
    get<3>(X).RealSpaceVectors = {{0., 0.}, {1., 1.}};

    get<0>(X).Variable = 2.;
    get<1>(X).Variable = 3;
    get<2>(X).Variable = -1;
    get<3>(X).LowerUpperBound = {0, 2};
    get<3>(X).Variable = 0;

    auto x = RealSpace(X);
    for(int j=0; j < 2; ++j){
      x[j] -= (j==0)*3-(j==1);
      assert(std::abs(x[j]) < 10*std::numeric_limits<Real>::epsilon() );
    }

    get<3>(X).Variable = 1;
    x = RealSpace(X);
    for(int j=0; j < 2; ++j){
      x[j] -= (j==0)*3-(j==1) + 1;
      assert(std::abs(x[j]) < 10*std::numeric_limits<Real>::epsilon() );
    }
    
  }

}//namespace
