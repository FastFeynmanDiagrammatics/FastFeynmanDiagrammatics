namespace ffd::relation_vector::unit_test{

  void generator(){
    std::vector<std::array<int, 2>> v;
    v.push_back({1, 2});
    v.push_back({1, 3});
    v.push_back({1, 1});
    v.push_back({4, 12});
    auto oracle = RelationOracle_g(v);
    using namespace ffd::user_space;
    std::cerr << oracle(3) << ' ' << oracle(4)
	      << oracle(1) << ' ' << oracle(20)
	      << oracle(2) << '\n';
  }

}//namespace
