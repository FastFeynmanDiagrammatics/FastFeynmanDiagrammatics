

namespace ffd::user_space::quadratic_action::unit_test{

  void creation(){
    using namespace std;

    
    int x1 = 1, x2 = 2;


    
    Real t = 1.;
    auto real_action_term = Bar(Psi_(1))(x1)*Psi_(1)(x2)*(-t);
    
    
    static_assert( is_same_v<decltype(real_action_term), QuadraticActionTerm<Real>> );
    assert( std::size(real_action_term) == 2 );
    assert( std::size(real_action_term[0]) == 1 );
    assert( std::size(real_action_term[1]) == 1 );
    assert( real_action_term[0][0] == Bar(Psi_(1))[0] );
    assert( std::any_cast<int>(real_action_term[0].Position) == x1 );
    assert( real_action_term[1][0] == Psi_(1)[0] );
    assert( std::any_cast<int>(real_action_term[1].Position) == x2 );
    assert( std::abs( real_action_term.Value + t ) < 10*std::numeric_limits<Real>::epsilon() );

    
    
    Complex z = 1.;
    auto complex_action_term = Bar(Psi_(1))(x1)*Psi_(1)(x2)*z;

    
    static_assert( is_same<decltype(complex_action_term), QuadraticActionTerm<Complex>>::value );
    assert( std::size(complex_action_term) == 2 );
    assert( std::size(complex_action_term[0]) == 1 );
    assert( std::size(complex_action_term[1]) == 1 );
    assert( complex_action_term[0][0] == Bar(Psi_(1))[0] );
    assert( std::any_cast<int>(complex_action_term[0].Position) == x1 );
    assert( complex_action_term[1][0] == Psi_(1)[0] );
    assert( std::any_cast<int>(complex_action_term[1].Position) == x2 );
    assert( std::abs( complex_action_term.Value - z ) < 10*std::numeric_limits<Real>::epsilon() );

    
    
    QuadraticAction real_action {Bar(Psi_(1))(x1)*Psi_(1)(x2)*(-t)};

    static_assert( is_same<decltype(real_action), QuadraticAction<Real>>::value );
    assert( std::size(real_action[0]) == 2 );
    assert( std::size(real_action[0][0]) == 1 );
    assert( std::size(real_action[0][1]) == 1 );
    assert( real_action[0][0][0] == Bar(Psi_(1))[0] );
    assert( std::any_cast<int>(real_action[0][0].Position) == x1 );
    assert( real_action[0][1][0] == Psi_(1)[0] );
    assert( std::any_cast<int>(real_action[0][1].Position) == x2 );
    assert( std::abs( real_action[0].Value + t ) < 10*std::numeric_limits<Real>::epsilon() );
    
    

    QuadraticAction complex_action (Bar(Psi_(1))(x1)*Psi_(1)(x2)*(z));

    static_assert( is_same<decltype(complex_action), QuadraticAction<Complex>>::value );
    
    
  }

}//namespace
