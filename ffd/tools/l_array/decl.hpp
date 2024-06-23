namespace ffd{

  template<typename T,
	   std::size_t n>
  struct l_array{
    std::vector<T> coef;

    l_array(): coef(n) {}

    template<typename T1>
    l_array(T1 x): coef(n, x) {}

    template<typename T1>
    void fill(T1 x){ for(std::size_t j=0; j<n; ++j) coef[j] = x;}

    inline T& operator[](std::size_t j){
      assert(( j>=0 && j<n ));
      return coef[j];
    }


    inline T const& operator[](std::size_t j) const{
      assert(( j>=0 && j<n ));
      return coef[j];
    }

  };

  
  template<typename T, std::size_t n>
  constexpr std::size_t
  size(l_array<T, n> const&){return n;}


  template<auto k, typename T, std::size_t n>
  T&
  get(l_array<T, n> & x){static_assert(k>=0 && k<n);
    return x[k];}


  template<auto k, typename T, std::size_t n>
  T const&
  get(l_array<T, n> const& x){static_assert(k>=0 && k<n);
    return x[k];}


  
  
}//namespace
