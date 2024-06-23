

namespace ffd::relation_sets::unit_test{

  template<typename T>
  void CheckPositionOfElement1InSetOfSets2(T const& element,
					   std::set<std::set<T>> const& set_of_sets,
					   int position){
    auto x = PositionOfElement1InSetOfSets2(element, set_of_sets);
    if(position<0){
      assert( !x.has_value() );
    }else{
      assert( x.has_value() );
      assert( x.value() == position);
    }
  }

  
  void position_of_element_1_in_set_of_sets_2(){

    std::set<std::set<char>> partition = {{0}, {1, 4}, {-4, -1}, {2, 3, -2}, {-3, 5, 8, 9}};
    assert(std::size(partition) == 5);
    
    CheckPositionOfElement1InSetOfSets2<char>(-4, partition, 0);
    CheckPositionOfElement1InSetOfSets2<char>(-1, partition, 0);
    CheckPositionOfElement1InSetOfSets2<char>(8, partition, 1);
    CheckPositionOfElement1InSetOfSets2<char>(1, partition, 4);
    CheckPositionOfElement1InSetOfSets2<char>(4, partition, 4);
    CheckPositionOfElement1InSetOfSets2<char>(-2, partition, 2);
    CheckPositionOfElement1InSetOfSets2<char>(0, partition, 3);
    CheckPositionOfElement1InSetOfSets2<char>(5, partition, 1);
    CheckPositionOfElement1InSetOfSets2<char>(10, partition, -1);
    CheckPositionOfElement1InSetOfSets2<char>(7, partition, -1);
    CheckPositionOfElement1InSetOfSets2<char>(-11, partition, -1);

  }

}//namespace
