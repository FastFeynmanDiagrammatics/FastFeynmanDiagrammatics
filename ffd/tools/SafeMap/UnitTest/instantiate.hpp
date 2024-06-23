namespace ffd::user_space::safe_map::unit_test{

  void instantiate(){

    SafeMap x(IntSequence<1, -1>(), std::array<char, 2>{'a', 'b'});
    auto y = CreateSafeMap<'a', 2>(std::array<std::string, 2>{"aa", "cc"});
    // std::cerr<<position_of_value_in_sequence<-1, 1, 2, 1>()<<std::endl;
    static_assert( static_position_of_value_in_sequence<1, 1, 2, 1>() == 0 );
    static_assert( static_position_of_value_in_sequence<-1, 1, 2, 1>() == -1 );
    [[maybe_unused]] auto yy = ffd::get<'a'>(y);
    // std::cerr<<ffd::get<-1>(y)<<std::endl;
    // std::cerr<<y[1]<<std::endl;
    // std::cerr<<ffd::get<0>(x)<<std::endl;

  }


}//namespace
