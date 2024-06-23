namespace ffd::periodic_coordinate::unit_test{

  void operator_parentheses_and_sign(){
    
    using ffd::sign;
    
    PeriodicCoordinate<int, 3> x;

    x.LowerUpperBound = {-1, 2};

    x.NotAntiPeriodic = false;

    x.Variable = 3;
    assert(x() == 0);
    assert(sign(x) == -1);

    x.Variable = 2;
    assert(x() == -1);
    assert(sign(x) == -1);
    
    x.Variable = -1;
    assert(x() == -1);
    assert(sign(x) == 1);

    x.Variable = -2;
    assert(x() == 1);
    assert(sign(x) == -1);

    x.Variable = -5;
    assert(x() == 1);
    assert(sign(x) == 1);

    
    PeriodicCoordinate<Real, 2> tau;

    tau.LowerUpperBound = {0., 1.};

    tau.NotAntiPeriodic = false;

    tau.Variable = 1.3;
    assert(std::abs(tau() - .3) < 10*std::numeric_limits<Real>::epsilon() );
    assert(sign(tau) == -1);

    
  }

}//namespace
