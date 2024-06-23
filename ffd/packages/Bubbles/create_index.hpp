

namespace ffd::bubbles{

  template<auto d>
  
  int
  
  create_index(std::array<int, d> x,
	       std::array<int, d> L){
    int index = 0, prod_L = 1;
    for( decltype(d) j=0; j<d; ++j){
      index += x[j]*prod_L;
      prod_L *= L[j];
    }

    return index;
  }
  

}//namespace
