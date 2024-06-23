namespace ffd{using uint = unsigned int;
  using ulong = unsigned long;
  namespace user_space {using uint = ffd::uint; using ulong = ffd::ulong;}}

namespace ffd::core_math{


  inline constexpr Real Pi = 3.1415926535897932384626433832795028841971693993l;

  
  Complex const I = Complex(Real(0), Real(1));

  
  constexpr Real Factorial(int n){
    return n == 0 ? 1 : n*Factorial(n-1);
  }


  template<typename IntType>
  constexpr IntType sqrt_int(IntType square){
    IntType square_root = 0;

    
    assert( square >= 0 );
    for(; square_root*square_root < square; ++square_root){}
    assert( square_root*square_root == square);
    return  square_root;
  }


  constexpr long
  pow_int(int base, int exponent){
    long ret = 1;
    for( int j=0; j<exponent; ++j){
      ret *= base;
    }
    return ret;
  }

  template <class int_t>
  inline auto
  popcount(int_t S)
  {
    return __builtin_popcount(S);
  }
  
}//namespace ffd::core_math





