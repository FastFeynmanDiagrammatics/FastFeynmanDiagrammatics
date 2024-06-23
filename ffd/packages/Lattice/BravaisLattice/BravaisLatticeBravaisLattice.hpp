namespace ffd::lattice{

  template<int d, char c>
  void InitializeBravais(BravaisLattice<d, c>&, Real){}
  
  template<int d, char c>
  void RotateXYBravais(BravaisLattice<d, c>&, Real){}

  
  template<int d, char c>
  BravaisLattice<d, c>::BravaisLattice(std::array<int, d> const& L, Real scale_, Real angle_in_pis_){
    for(int j=0; j<d; ++j){
      auto& [low, upp] = LowerUpperBounds[j];
      low = -((L.at(j)+1)/2)+1;
      upp = L.at(j)/2+1;
    }
    InitializeBravais(*this, scale_);
    RotateXYBravais(*this, angle_in_pis_);
  }

  
}//namespace
