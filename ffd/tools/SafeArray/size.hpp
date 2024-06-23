namespace ffd::user_space{

  template<typename value_t, std::size_t size_v>
  inline constexpr std::size_t size(SafeArray<value_t, size_v> const&){return size_v;}

}//namespace
