

namespace ffd::wick_matrix::unit_test{

  using namespace ffd::user_space;

  using namespace ffd::qf_product;
  
  void Creation_qf_product(QuantumFieldProduct P_){
    using std::size;
    WickMatrix<Real> M0(P_[0], 0);
    WickMatrix<Real> M0P(P_, 0);
    assert(size(M0) == 0);
    assert(size(M0P) == 0);
    WickMatrix<Real> M1(P_[0], 1);
    WickMatrix<Real> M1P(P_, 1);
    assert(size(M1) == 1);
    assert(size(M1P) == 1);
    WickMatrix<Real> M10(P_[0], 10);
    WickMatrix<Real> M10P(P_, 10);
    assert(size(M10) == 10);
    assert(size(M10P) == 10);
    if(P_[0].Dagger() != 0){
      assert(size(M10.Components) == 100);
    }else if(P_[0].IsFermion() && P_[0].Dagger() == 0){
      assert(size(M10.Components) == 45);
    }else{
      assert(size(M10.Components) == 55);
    }
  }


  void Creation(){
    Creation_qf_product(Psi_(1));
    Creation_qf_product(Phi_(1));
    Creation_qf_product(Rho_(1));
    Creation_qf_product(Eta_(1));
  }

  
}//namespace
