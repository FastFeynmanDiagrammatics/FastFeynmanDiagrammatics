

namespace ffd::qf_vertex::unit_test{

  void Creation(){
    QuantumFieldDot D_;
    QuantumFieldVertex V_(D_);

    assert(V_[0] == D_);

    QuantumFieldVertex V2_;

    assert(std::size(V2_) == 0);

    V2_.push_back(D_);

    assert(V2_[0] == D_);
  }
  
}
