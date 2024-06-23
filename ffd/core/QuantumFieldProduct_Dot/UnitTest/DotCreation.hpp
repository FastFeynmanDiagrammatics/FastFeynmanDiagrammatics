

namespace ffd::qf_dot::unit_test{

  void DotCreation(){

    QuantumFieldDot D_;
    D_.resize(10);
    D_.Position = static_cast<float>(1);
    
    assert(std::size(D_) == 10);
  }

}
