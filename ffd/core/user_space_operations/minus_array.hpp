namespace ffd::user_space{

  template<typename T, auto d>
  auto
  operator-(std::array<T, d> const& x){
    auto y = x;
    for(uint j=0; j<d; ++j){
      y[j] = - y[j];
    }
    return y;
  }

}//namespace
