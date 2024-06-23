namespace ffd::periodic_coordinate{

  void periodic_coordinates(){
    
    using ffd::sign, ffd::get;
    
    PeriodicCoordinates<3, int, double> X;

    get<0>(X).NotAntiPeriodic = false;
    get<0>(X).LowerUpperBound = {0, 2};
    get<0>(X).Variable = 2;

    assert( sign(get<0>(X)) == -1);
    assert( sign(get<1>(X)) ==  1);
    assert( sign(X)         == -1);
    
  }

}//namespace
