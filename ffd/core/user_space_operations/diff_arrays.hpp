namespace ffd::user_space{

  template<typename T, auto d>
  auto
  operator-(std::array<T, d> const& x0,
	    std::array<T, d> const& x1){
    auto x2 = x0;
    for(uint j=0; j<d; ++j){
      x2[j] -= x1[j];
    }
    return x2;

  }

}//namespace
