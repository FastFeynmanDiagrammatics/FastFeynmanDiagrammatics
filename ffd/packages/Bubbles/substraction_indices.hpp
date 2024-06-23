

namespace ffd::bubbles{

  template<int d, typename T1, typename T2, typename T3>
  
  std::array<int, d>
  
  substraction_indices(T1 a,
		       T2 b,
		       T3 L){
    using ffd::vector_range::Range;
    
    auto a_b = a;
    for( int j: Range(d) ){
      a_b[j] -= b[j];
      a_b[j] += a_b[j] < 0 ? L[j] : 0;
    }

    return a_b;
  }


}//namespace
