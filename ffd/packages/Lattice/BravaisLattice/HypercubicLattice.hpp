namespace ffd::lattice{

  template<int d>
  using HypercubicLattice = BravaisLattice<d, 'H'>;

  template<int d>
  void InitializeBravais(HypercubicLattice<d>& L_, Real LatticeConstant){
    for(int j=0; j < d; ++j){
      auto& BravaisVector = L_.BravaisVectors[j];
      for(int k=0; k < d; ++k){
	BravaisVector[k] = LatticeConstant*(j==k);
      }
    }
  }
  
}//namespace ffd::lattice
