namespace ffd::lattice{

  template<int space_dimensions, int num_atoms_unit_cell, char unit_cell_name = '0'>
  class UnitCell{
  public:
        
    std::array<std::array<Real, space_dimensions>, num_atoms_unit_cell> AtomicOffsets;

    UnitCell(Real scale_factor = 1.);
    
  };

}//namespace ffd::lattice
