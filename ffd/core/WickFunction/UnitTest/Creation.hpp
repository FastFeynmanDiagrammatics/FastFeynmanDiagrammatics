

namespace ffd::wick_function::unit_test{

  using namespace ffd::user_space;
  
  void Creation(){
    WickFunction W;
    assert(W.WhichBlock(Rho_(0)[0]) == 31);
    assert(W.WhichBlock(Rho_(1)[0]) == 32);
  }
  
}//namespace
