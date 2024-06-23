

namespace ffd::relation_sets::unit_test{

  void merge(){

    std::set<char> S1 = {1, 3, 5, 3, 4, 7};
    assert(std::size(S1) == 5);

    std::set<char> S2 = {1, 2, 4, 3, 8, 1, 3};
    assert(std::size(S2) == 5);

    auto S3 = Merge(S1, S2);
    assert(std::size(S3) == 7);
    std::set<char> S3_eq = {1, 2, 3, 4, 5, 7, 8};
    assert(S3 == S3_eq);

    
  }

}//namespace
