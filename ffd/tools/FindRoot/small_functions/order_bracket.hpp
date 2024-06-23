

namespace ffd::find_root{

  std::array<Real, 2>

  order_bracket(std::array<Real, 2> x){
    if( x[0] > x[1] ){
      std::swap(x[0], x[1]);
    }
    
    return x;
  }


}//namespace
