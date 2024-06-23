

namespace ffd::relation_sets::unit_test{

  void creation(){

    RelationSets<int> R;
    std::set<int> S = {1, 2, 3};
    R.insert(S);
    R = {{1, 2, 3}, {3}, {5, 6}, {}};

    Element0IsInRelationWithElement1<int> x_tilde_y;
    x_tilde_y = {3, 5};
    
  }

}//namespace
