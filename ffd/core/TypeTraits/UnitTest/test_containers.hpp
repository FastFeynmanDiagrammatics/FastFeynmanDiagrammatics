namespace ffd::type_traits::unit_test{

void test_containers(){
    static_assert(!(is_container<int>::value));
    static_assert(!(is_container<double>::value));
    static_assert(!(is_container<Complex>::value));
    static_assert(is_container<std::array<double,1> >::value);
    static_assert(is_container<std::array<double,3> >::value);
    static_assert(is_container<std::vector<double> >::value);
    static_assert(is_container<std::vector<double> >::value);

}

}//namespace
