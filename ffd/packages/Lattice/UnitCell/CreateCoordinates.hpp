

namespace ffd::lattice{

  template<int d, int n_atoms, char c>
  auto CreateCoordinates(UnitCell<d, n_atoms, c> const& U_){
    using ffd::get;
    ffd::periodic_coordinate::PeriodicCoordinates<d, int> ret;
    get<0>(ret).LowerUpperBound = {0, n_atoms};
    get<0>(ret).RealSpaceVectors.resize(n_atoms);
    for(int j=0; j < n_atoms; ++j){
      get<0>(ret).RealSpaceVectors[j] = U_.AtomicOffsets[j];
    }
    get<0>(ret).Variable = 0;
    return ret;
  }

}//namespace ffd::lattice
