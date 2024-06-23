namespace ffd::relation_vector::unit_test{

  void relation_test(){
    Container<int> C;
    AddRelation_r({1, 2},
		  C);
    AddRelation_r({3, 3},
		  C);
    AddRelation_r({3, 3},
		  C);
    AddRelation_r({2, 1},
		  C);
    AddRelation_r({4, -20},
		  C);
    AddRelation_r({4, 3},
		  C);
    AddRelation_r({2, 2},
		  C);
    
    std::cerr << C << '\n';
  }

}//namespace
