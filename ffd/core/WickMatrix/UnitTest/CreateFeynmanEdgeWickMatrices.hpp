

namespace ffd::wick_matrix::unit_test{

  using namespace ffd::quantum_field;
  
  template<typename T, typename intType = intTypeDefault>
  void CheckElementsWickMatrix(WickMatrix<T, intType> const& M,
			       std::array<char, 2> index,
			       ffd::quantum_field::QuantumField QF0_,
			       std::array<char, 2> VD0_,
			       ffd::quantum_field::QuantumField QF1_,
			       std::array<char, 2> VD1_){
    auto [QFVD0, QFVD1] = M.Components[ index[1] + index[0]*size(M) ];
    auto [QF0, VD0] = QFVD0;
    auto [QF1, VD1] = QFVD1;
    assert(QF0 == QF0_);
    assert(QF1 == QF1_);
    assert(VD0 == VD0_);
    assert(VD1 == VD1_);
  }

  
  void CreateFeynmanEdgeWickMatrices(){
    using namespace ffd::user_space;
    using std::size;
    int x = 1, y=2, z=3;
    auto D = (Psi_(1)*Psi_(-1)*Bar(Eta_(0)))(x);
    auto V = Bar(D)*D;
    auto Green = Psi_(1)(y)*Bar(Psi_(1))(z);
    auto G = V | V | V | Green;
    WickFunction W;
    for(std::size_t j=0; j < size(G); ++j){ 
      W *= G(j);
    }
    
    assert(W.Sign == -1);
    
    auto Matrices = ffd::user_space::CreateFeynmanEdgeWickMatrices<int>(W);

    for(char j=0; j < char(size(Matrices[0])); ++j){
      for(char k=0; k < char(size(Matrices[0])); ++k){
	CheckElementsWickMatrix(Matrices[0], {j, k}, Bar(Eta_(0)[0]), {j, 1}, Eta_(0)[0], {k, 0});
      }
    }

    for(char j=0; j < char(size(Matrices[1])); ++j){
      for(char k=0; k < char(size(Matrices[1])); ++k){
	CheckElementsWickMatrix(Matrices[1], {j, k}, Bar(Psi_(-1)[0]), {j, 0}, Psi_(-1)[0], {k, 1});
      }
    }


    for(char j=0; j < char(size(Matrices[1])); ++j){
      for(char k=0; k < char(size(Matrices[1])); ++k){
	CheckElementsWickMatrix(Matrices[2], {j, k}, Bar(Psi_(1)[0]), {j, 0}, Psi_(1)[0], {k, 1});
      }
    }

    for(char j=0; j < char(size(Matrices[1])); ++j){
      char k = static_cast<char>(char(size(Matrices[1])));
      CheckElementsWickMatrix(Matrices[2], {j, k}, Bar(Psi_(1)[0]), {j, 0}, Psi_(1)[0], {k, 0});
    }

    char j = static_cast<char>(size(Matrices[1]));
    for(char k=0; k < char(size(Matrices[1])); ++k){
      CheckElementsWickMatrix(Matrices[2], {j, k}, Bar(Psi_(1)[0]), {j, 1}, Psi_(1)[0], {k, 1});
    }
    
    CheckElementsWickMatrix(Matrices[2], {j, j}, Bar(Psi_(1)[0]), {j, 1}, Psi_(1)[0], {j, 0});
    
  }

  
}//namespace
