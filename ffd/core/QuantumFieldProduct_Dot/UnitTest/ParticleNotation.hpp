


namespace ffd::qf_product::unit_test{

  
  using namespace ffd::user_space;

  
  void ParticleNotation(){
    auto P_up = Psi_("up");
    auto P_down = Psi_("down");
    auto P1 = Psi_(1);
    
    static_assert(std::is_same<decltype(P_up), decltype(P1)>::value);
    assert(P_up == P1);
    assert(P_up == FlipSpin(P_down));
    assert(P1 != P_down);

    
    auto P2 = Phi_(1);
    auto P3 = Eta_(1);
    auto P4 = Rho_(1);
    auto P5 = Nambu(Psi_(1));
    auto P6 = Nambu(Psi_(-1));
    

    assert(P4 == Bar(P4));
    assert(P2 == Bar(P2));
    assert(Psi_(0) == FlipSpin(Psi_(0)));
    assert(Bar(P1*P3) == Bar(P3)*Bar(P1));
    assert(FlipSpin(P1*P3) == FlipSpin(P1)*FlipSpin(P3));
    assert(P5 != P1);
    assert(FlipSpin(P1) != P6);
    assert(!P5[0].NotNambu());
    assert(!P6[0].NotNambu());
    
    static_assert(std::is_same<decltype(P1), decltype(P3)>::value);
    static_assert(std::is_same<decltype(P1), decltype(P4*P2)>::value);

    
    double x=1;
    static_assert(std::is_same<decltype(P1(x)), ffd::qf_dot::QuantumFieldDot>::value);
    static_assert(std::is_same<decltype((P1*P2*P3)(x)), ffd::qf_dot::QuantumFieldDot>::value);
    static_assert(std::is_same<decltype((FlipSpin(P1*P2)*Bar(P3*P4))(x)),ffd::qf_dot::QuantumFieldDot>::value);

    
    assert(P1 != P2);
    assert(P3 != P4);
    assert(P3 < P1);
    assert(P4 < P1);
    assert(P2 < P4);

    
    assert(P3*P4 < P1*P1);
    assert(P3*P4 < P1*P1*P1);
    
    auto P10 = P1*P2*P3;

    assert(P10[0] == P1[0]);
    assert(P10[1] == P2[0]);
    assert(P10[2] == P3[0]);
    
  }

}//namespace ffd::qf_product::unit_test
