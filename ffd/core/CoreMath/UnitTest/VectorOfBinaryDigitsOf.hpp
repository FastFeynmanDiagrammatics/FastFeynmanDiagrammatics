

namespace ffd::set_theory::unit_test{

  void CheckVectorBinaryDigits(){
    auto v0 = VectorOfBinaryDigitsOf(0);
    assert(std::size(v0) == 0);

    auto v1 = VectorOfBinaryDigitsOf(0b10);
    assert(std::size(v1) == 1);
    assert(v1[0] == 1);

    auto v2 = VectorOfBinaryDigitsOf(0b101);
    assert(std::size(v2) == 2);
    for(std::size_t j=0; j < std::size(v2); ++j){
      assert(v2[j] == 0*(j==0)+2*(j==1));
    }

    auto v3 = VectorOfBinaryDigitsOf(0b110101);
    assert(std::size(v3) == 4);
    for(std::size_t j=0; j < std::size(v3); ++j){
      assert(v3[j] == 0*(j==0)+2*(j==1)+4*(j==2)+5*(j==3));
    }

    auto v4 = VectorOfBinaryDigitsOf(0b1101010);
    assert(std::size(v4) == 4);
    for(std::size_t j=0; j < std::size(v4); ++j){
      assert(v4[j] == 1*(j==0)+3*(j==1)+5*(j==2)+6*(j==3));
    }
    
  }

}
