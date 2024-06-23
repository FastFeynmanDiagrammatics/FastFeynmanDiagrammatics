namespace ffd::user_space::operations::unit_test{

  void sum(){
    std::tuple<int> t0 , t1;
    std::get<0>(t0) = 2;
    std::get<0>(t1) = 3;
    auto t2 = t0+t1;
    std::cerr << t2 << '\n';
    std::tuple<std::array<int, 2>, int> x0, x1;
    std::get<0>(x0) = {1, 2};
    std::get<0>(x1) = {-1, 8};
    std::get<1>(x0) = 3;
    std::get<1>(x1) = 7;
    auto const x2 = x0 + x1;
    std::cerr << x2;
  }

}//namespace
