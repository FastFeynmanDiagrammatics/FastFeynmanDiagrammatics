namespace ffd::lattice{

  template<int d>
  using SegmentUnitCell = UnitCell<d, 2, 'S'>;

  template<int d>
  void InitializeUnitCell(SegmentUnitCell<d>& S_, Real Scale){
    S_.AtomicOffsets[0] = {0.};
    S_.AtomicOffsets[1] = {0.};
    S_.AtomicOffsets[1][0] = Scale;
  }

}
