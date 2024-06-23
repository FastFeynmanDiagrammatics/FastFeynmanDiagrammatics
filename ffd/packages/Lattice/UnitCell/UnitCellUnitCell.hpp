

namespace ffd::lattice{

  template<int j, int k, char c>
  void InitializeUnitCell(UnitCell<j, k, c>&, Real) {}

  template<int d, int n_atoms, char name>
  UnitCell<d, n_atoms, name>::UnitCell(Real Scale_){
    InitializeUnitCell(*this, Scale_);
  }

}//namespace
