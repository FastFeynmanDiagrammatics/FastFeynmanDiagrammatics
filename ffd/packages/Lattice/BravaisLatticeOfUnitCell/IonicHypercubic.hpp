namespace ffd::lattice{

  template<int d>
  auto
  IonicHypercubic(std::array<int, d> const& L){
    BravaisLatticeOfUnitCell<d, 'H', 2, 'S'> ret(L);
    for(uint j=0; j<d; ++j){
      ret.BravaisLatticeObject.BravaisVectors[j][0] += 1.;
    }
    return ret;
  }
  
}//namespace ffd::lattice
