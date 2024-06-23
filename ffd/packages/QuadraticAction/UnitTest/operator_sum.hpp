

namespace ffd::user_space::quadratic_action::unit_test{

  void operator_sum(){
    
    using namespace std;


    std::array<int, 2> x0{0, 0}, x1{1, 0};
    Real t = 1;
    QuadraticAction action(Bar(Psi_(1))(x1)*Psi_(1)(x0)*(-t));
    action += Bar(Psi_("down"))(x1)*Psi_("down")(x0)*(-t);
    
    action += HermitianConjugate(action);


    
    Complex t2{1., 1.};
    QuadraticAction action_c( Bar(Psi_("up"))(x1)*Psi_("up")(x0)*t2 );
    action_c += Bar(Psi_("down"))(x1)*Psi_("down")(x0)*t2;

    action_c += HermitianConjugate(action_c);

    assert( abs( imag( action_c[0].Value + action_c[3].Value ) ) < numeric_limits<Real>::epsilon() );
    assert( abs( imag( action_c[1].Value + action_c[2].Value ) ) < numeric_limits<Real>::epsilon() );
    assert( abs( real( action_c[1].Value - action_c[2].Value ) ) < numeric_limits<Real>::epsilon() );

    auto x3 = any_cast<array<int, 2>>( action_c[0][0].Position );
    assert( x3  == x1 );
    auto x4 = any_cast<array<int, 2>>( action_c[0][1].Position );
    assert( x4  == x0 );
    auto x5 = any_cast<array<int, 2>>( action_c[3][1].Position);
    assert( x5  == x1 );
    auto x6 = any_cast<array<int, 2>>( action_c[2][0].Position);
    assert( x6  == x0 );
    assert( action_c[3][0] == Bar(Psi_("down")) );
    assert( action_c[2][1] == Psi_("up") );
    
  }

}//namespace
