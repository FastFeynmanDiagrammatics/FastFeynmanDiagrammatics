

namespace ffd::wick_function::unit_test{

  using namespace ffd::user_space;
  
  void operator_multiply(){
    using std::size;
    using namespace ffd::quantum_field;
    using namespace ffd::feynman_edge;
    int x = 123;
    auto D = (Psi_(-1)*Psi_(1)*Bar(Eta_(0)))(x);
    auto BD = Bar(D);
    auto DHubbard = Psi_(-1)*Psi_(1);
    auto H = (Bar(DHubbard)*DHubbard)(x);
    auto G = H | H | BD*D | BD*D*BD*D;

    WickFunction W;
    W *= G(0);

    assert(size(W) == 2);
    for(auto j: {0, 1}){
      for(auto k: {0, 1}){
	assert(size(W[j][k]) == 1);
      }
    }
    auto Psi_up = Psi_(1)[0];
    auto Psi_down = Psi_(-1)[0];
    std::array<char, 2> VertexDot = {0, 0};
    std::pair<QuantumField, std::array<char, 2>> QFVD_barup = {Bar(Psi_up), VertexDot};
    std::pair<QuantumField, std::array<char, 2>> QFVD_bardown = {Bar(Psi_down), VertexDot};
    std::pair<QuantumField, std::array<char, 2>> QFVD_up = {Psi_up, VertexDot};
    std::pair<QuantumField, std::array<char, 2>> QFVD_down = {Psi_down, VertexDot};
    auto WQFVD = W[0][0][0];
    static_assert(std::is_same<QuantumFieldVertexDot, decltype(WQFVD)>::value);
    assert(WQFVD == QFVD_barup);
    WQFVD = W[0][1][0];
    assert(WQFVD == QFVD_up);
    WQFVD = W[1][0][0];
    assert(WQFVD == QFVD_bardown);
    WQFVD = W[1][1][0];
    assert(WQFVD == QFVD_down);
    assert(size(W.BlockNumbers) == 2);
    assert(size(W.IsFermionicBlock) == 2);
    for(auto j: {0, 1}){
      assert(W.BlockNumbers[j] == W.WhichBlock(Psi_(1-2*j)[0]));
      assert(W.IsFermionicBlock[j]);
    }
    assert(W.Sign == 1);

    
    W *= G(1);

    assert(size(W) == 2);
    for(auto j: {0, 1}){
      for(auto k: {0, 1}){
	assert(size(W[j][k]) == 2);
      }
    }
    assert(size(W.BlockNumbers) == 2);
    assert(size(W.IsFermionicBlock) == 2);
    for(auto j: {0, 1}){
      assert(W.BlockNumbers[j] == W.WhichBlock(Psi_(1-2*j)[0]));
      assert(W.IsFermionicBlock[j]);
    }
    assert(W.Sign == 1);


    W *= G(2);

    assert(size(W) == 3);
    for(auto j: {0, 1}){
      for(auto k: {0, 1}){
	assert(size(W[j][k]) == 3);
      }
    }
    for(auto j: {0, 1}){
      assert(size(W[2][j]) == 1);
    }
    assert(size(W.BlockNumbers) == 3);
    assert(size(W.IsFermionicBlock) == 3);
    for(auto j: {0, 1}){
      assert(W.BlockNumbers[j] == W.WhichBlock(Psi_(1-2*j)[0]));
      assert(W.IsFermionicBlock[j]);
    }
    assert(W.BlockNumbers[2] == W.WhichBlock(Eta_(0)[0]));
    assert(!W.IsFermionicBlock[2]);
    assert(W.Sign == 1);


    W *= G(3);

    assert(size(W) == 3);
    for(auto j: {0, 1}){
      for(auto k: {0, 1}){
	assert(size(W[j][k]) == 5);
      }
    }
    for(auto j: {0, 1}){
      assert(size(W[2][j]) == 3);
    }
    QFVD_barup.second = {3, 0};
    WQFVD = W[0][0][3];
    assert(WQFVD == QFVD_barup);
    QFVD_barup.second = {3, 2};
    WQFVD = W[0][0][4];
    assert(WQFVD == QFVD_barup);
    QFVD_up.second = {3, 1};
    WQFVD = W[0][1][3];
    assert(WQFVD == QFVD_up);
    QFVD_up.second = {3, 3};
    WQFVD = W[0][1][4];
    assert(WQFVD == QFVD_up);

    QFVD_bardown.second = {3, 0};
    WQFVD = W[1][0][3];
    assert(WQFVD == QFVD_bardown);
    QFVD_bardown.second = {3, 2};
    WQFVD = W[1][0][4];
    assert(WQFVD == QFVD_bardown);
    QFVD_down.second = {3, 1};
    WQFVD = W[1][1][3];
    assert(WQFVD == QFVD_down);
    QFVD_down.second = {3, 3};
    WQFVD = W[1][1][4];
    assert(WQFVD == QFVD_down);

    auto Eta = Eta_(0)[0];
    std::pair<QuantumField, std::array<char, 2>> QFVD_eta    = {Eta,      QFVD_up.second};
    std::pair<QuantumField, std::array<char, 2>> QFVD_bareta = {Bar(Eta), QFVD_up.second};
    QFVD_bareta.second = {3, 1};
    WQFVD = W[2][0][1];
    assert(WQFVD == QFVD_bareta);
    QFVD_bareta.second = {3, 3};
    WQFVD = W[2][0][2];
    assert(WQFVD == QFVD_bareta);
    QFVD_eta.second = {3, 0};
    WQFVD = W[2][1][1];
    assert(WQFVD == QFVD_eta);
    QFVD_eta.second = {3, 2};
    WQFVD = W[2][1][2];
    assert(WQFVD == QFVD_eta);

    assert(size(W.BlockNumbers) == 3);
    assert(size(W.IsFermionicBlock) == 3);
    for(auto j: {0, 1}){
      assert(W.BlockNumbers[j] == W.WhichBlock(Psi_(1-2*j)[0]));
      assert(W.IsFermionicBlock[j]);
    }
    assert(W.BlockNumbers[2] == W.WhichBlock(Eta_(0)[0]));
    assert(!W.IsFermionicBlock[2]);
    assert(W.Sign == 1);

    
  }




}//namespace
