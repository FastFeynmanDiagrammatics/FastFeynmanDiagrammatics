namespace ffd::user_space{

  template<typename T>
  T Nambu(T const& P_){
    T ret = P_;
    for(std::size_t j=0; j < size(P_); ++j){
      ret[j] = Nambu( ret[j] );
    }
    return ret;
  }

}//namespace
