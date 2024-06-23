

namespace ffd::flat_map::unit_test{

  void UnitTest(){

    FlatMap<int, int> Map;

    Map.insert_or_assign(2, 3);
    Map.insert_or_assign(1, 2);
    Map.insert_or_assign(6, 7);
    Map.insert_or_assign(9, 10);
    Map.insert_or_assign(3, 4);
    Map.insert_or_assign(-1, 0);
    
    assert(Map.at(6) == 7);
    assert(Map.at(2) == 3);
    assert(Map.at(1) == 2);

    assert(Map.count(6) == 1);
    assert(Map.count(0) == 0);
    assert(Map.count(-1) == 1);

    Map.insert_or_assign(2, -3);

    assert(Map.at(2) == -3);

    
  }
  
}
