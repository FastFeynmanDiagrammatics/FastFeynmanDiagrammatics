

namespace ffd::wick_block_oracle::unit_test{

  void UnitTest(){
    using namespace ffd::user_space;
    
    int x1=1, x2= 2;
    Real t = 1.;
    QuadraticAction<Real> A(Bar(Psi_("up"))(x1)*Psi_("up")(x2)*t);
    A += FlipSpin(A);

    WickBlockOracle oracle(A);

    assert( oracle( Psi_("down")[0] ) == 0 );
    assert( oracle( Bar(Psi_("up")[0]) ) == 1 );


    
    QuadraticAction<Complex> A2(Bar(Eta_(1))(x1)*Eta_(2)(x2)*Complex{1, 1});
    A2 += FlipSpin(A2);
    A2 += Bar(Eta_(1))(x1)*Eta_(3)(x2)*Complex{-1, -1};

    WickBlockOracle oracle2(A2);

    assert( oracle2( Eta_(-1)[0] ) == 0 );
    assert( oracle2( Eta_(-2)[0] ) == 0 );
    assert( oracle2( Bar(Eta_(2)[0]) ) == 1 );
    assert( oracle2( Eta_(3)[0] ) == 1 );
    

  }

}//namespace
