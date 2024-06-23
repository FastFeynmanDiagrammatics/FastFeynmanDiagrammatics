namespace ffd::type_traits{

template <typename T>
struct is_container {
  static const bool value = false;
};

template <typename T,typename Alloc>
struct is_container<std::vector<T,Alloc> > {
  static const bool value = true;
};

template <typename T, std::size_t S>
struct is_container<std::array<T,S> > {
  static const bool value = true;
};

}//namespace
