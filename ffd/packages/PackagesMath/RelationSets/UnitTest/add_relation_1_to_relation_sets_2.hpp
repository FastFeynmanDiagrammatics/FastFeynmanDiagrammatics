namespace ffd::relation_sets::unit_test{
  
  template<typename T>
  struct CheckAddRelation1ToRelationSets2{
    
    RelationSets<T> relation_sets;

    void AddRelation(std::array<T, 2> const& relation_elements,
		     RelationSets<T> const& relation_sets_after_){
      AddRelation1ToRelationSets2<T>(relation_elements, relation_sets);
      assert(relation_sets == relation_sets_after_);
    }
    
  };

  
  void add_relation_1_to_relation_sets_2(){

    CheckAddRelation1ToRelationSets2<char> R;
    
    R.AddRelation({1, 1}, {{1}});

    R.AddRelation({-1, -1}, {{1}, {-1}});

    R.AddRelation({2, 3}, {{1}, {-1}, {2, 3}});

    R.AddRelation({0, 1}, {{0, 1}, {-1}, {2, 3}});

    R.AddRelation({0, 0}, {{0, 1}, {-1}, {2, 3}});

    R.AddRelation({0, 1}, {{0, 1}, {-1}, {2, 3}});

    R.AddRelation({0, -1}, {{0, 1, -1}, {2, 3}});

    R.AddRelation({3, 4}, {{0, 1, -1}, {2, 3, 4}});

    R.AddRelation({1, -1}, {{0, 1, -1}, {2, 3, 4}});

    R.AddRelation({2, 5}, {{0, 1, -1}, {2, 3, 4, 5}});

    R.AddRelation({-3, -5}, {{-3, -5}, {0, 1, -1}, {2, 3, 4, 5}});

    R.AddRelation({1, 4}, {{-3, -5}, {0, 1, -1, 2, 3, 4, 5}});

    R.AddRelation({-3, -7}, {{-3, -5, -7}, {0, 1, -1, 2, 3, 4, 5}});

    R.AddRelation({10, 5}, {{-3, -5, -7}, {0, 1, -1, 2, 3, 4, 5, 10}});

    R.AddRelation({-2, 15}, {{-3, -5, -7}, {-2, 15}, {0, 1, -1, 2, 3, 4, 5, 10}});

    R.AddRelation({-10, -15}, {{-10, -15}, {-3, -5, -7}, {-2, 15}, {0, 1, -1, 2, 3, 4, 5, 10}});

    R.AddRelation({-100, -90}, {{-100, -90}, {-10, -15}, {-3, -5, -7},
				{-2, 15}, {0, 1, -1, 2, 3, 4, 5, 10}});

    R.AddRelation({-100, -90}, {{-100, -90}, {-10, -15}, {-3, -5, -7},
				{-2, 15}, {0, 1, -1, 2, 3, 4, 5, 10}});

    R.AddRelation({-10, -90}, {{-100, -90, -10, -15}, {-3, -5, -7},
			       {-2, 15}, {0, 1, -1, 2, 3, 4, 5, 10}});

    R.AddRelation({-7, 3}, {{-100, -90, -10, -15}, {-2, 15}, {-3, -5, -7, 0, 1, -1, 2, 3, 4, 5, 10}});

    R.AddRelation({-7,-15}, {{-2, 15}, {-100, -90, -10, -15, -3, -5, -7, 0, 1, -1, 2, 3, 4, 5, 10}});

    R.AddRelation({-2,-15}, {{-2, 15, -100, -90, -10, -15, -3, -5, -7, 0, 1, -1, 2, 3, 4, 5, 10}});
  }
  

}//namespace
