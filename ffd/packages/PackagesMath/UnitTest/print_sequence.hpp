

namespace ffd::packages_math::unit_test{

    
  void print_sequence(){
    StaticRangeSequence<1, 6> seq_1_6;
    auto vec = PrintSequenceRange(seq_1_6);
    assert(std::size(vec) == 5);
    for(int j=0; j < 5; ++j){
      assert(vec[j] == j+1);
    }
  }

}//namespace
