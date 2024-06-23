

namespace ffd::relation_sets::unit_test{

  void position_of_element_1_in_set_2(){

    std::set<char> S = {1, 2, 12, -4, 22, 3, 8};

    assert( !PositionOfElement1InSet2<char>(4, S).has_value() );
    assert( !PositionOfElement1InSet2<char>(-9, S).has_value() );
    
    assert( PositionOfElement1InSet2<char>(2, S).has_value() );
    assert( PositionOfElement1InSet2<char>(2, S).value() == 2);

    assert( PositionOfElement1InSet2<char>(-4, S).has_value() );
    assert( PositionOfElement1InSet2<char>(-4, S).value() == 0);

    
  }


}//namespace
