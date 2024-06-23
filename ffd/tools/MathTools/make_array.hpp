namespace ffd{

  template<typename T1, typename T2>
  auto
  make_array(T1 x1, T2 x2){
    static_assert(std::is_same_v<decltype(x1), decltype(x2)>);

    
    return std::array<decltype(x1), 2>{{x1, x2}};
  }

	
}//namespace
