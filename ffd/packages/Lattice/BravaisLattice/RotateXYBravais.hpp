

namespace ffd::lattice{

  template<char c>
  void RotateXYBravais(BravaisLattice<2, c>& L_, Real angle_in_pis_){
    std::array<std::array<Real, 2>, 2> V_Mat;
    for(auto j: {0, 1}){
      for(auto k: {0, 1}){
	V_Mat[j][k] = L_.BravaisVectors[j][k];
      }
    }
    Real Cos_a = cos(angle_in_pis_*ffd::core_math::Pi);
    Real Sin_a = sin(angle_in_pis_*ffd::core_math::Pi);
    for(auto j: {0, 1}){
      L_.BravaisVectors[j][0] = Cos_a*V_Mat[j][0] - Sin_a*V_Mat[j][1];
      L_.BravaisVectors[j][1] = Sin_a*V_Mat[j][0] + Cos_a*V_Mat[j][1];
    }
  }

}//namespace ffd::lattice
