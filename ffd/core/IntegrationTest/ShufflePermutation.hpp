

namespace ffd::core_integration_test{

  std::random_device rd;  //Will be used to obtain a seed for the random number engine                            
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

  void ShufflePermutation(std::vector<char>& permutation_){
    using std::size;
    std::uniform_int_distribution<char> dis(0, size(permutation_)-1);
    char j1 = dis(gen), j2 = dis(gen);
    while(j1 == j2){
      j2 = dis(gen);
    }
    char temp = permutation_[j1];
    permutation_[j1] = permutation_[j2];
    permutation_[j2] = temp;
  }

}//namespace
