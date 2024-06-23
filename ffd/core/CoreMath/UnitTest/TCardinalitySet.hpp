

namespace ffd::set_theory::unit_test{

  void TestCardinalitySet(){

    assert( CardinalitySet(1) == 1);
    assert( CardinalitySet(0) == 0);
    assert( CardinalitySet(3) == 2);
    assert( CardinalitySet(0b010011011010) == 6);
    assert( CardinalitySet(0b0101110110101) == 8);
    assert( CardinalitySet(0b0111010111101) == 9);
    
  }

}
