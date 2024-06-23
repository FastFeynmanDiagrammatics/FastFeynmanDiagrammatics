

namespace ffd::lattice{

  template<int d>
  using EquilateralTriangleUnitCell = UnitCell<d, 3, 'E'>;

  template<int d>
  void InitializeUnitCell(EquilateralTriangleUnitCell<d>& E_, Real Scale){
    E_.AtomicOffsets[0] = {0.};
    E_.AtomicOffsets[1] = {0.};
    E_.AtomicOffsets[1][0] = Scale;
    E_.AtomicOffsets[2] = {0.};
    E_.AtomicOffsets[2][0] = 0.5l*Scale;
    E_.AtomicOffsets[2][1] = 0.86602540378443864679l*Scale;
  }
  
}//namespace ffd::lattice
