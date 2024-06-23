namespace ffd::lattice{
  
  template<int d, char bravais_name, int atoms_unit_cell, char unit_cell_name>
  class BravaisLatticeOfUnitCell{
  public:

    static constexpr int dimension = d;

    static constexpr int number_atoms_unit_cell = atoms_unit_cell;
    

    static Real ScaleBravais;

    static Real AngleInGreekPis;

    
    BravaisLattice<d, bravais_name> BravaisLatticeObject;

    UnitCell<d, atoms_unit_cell, unit_cell_name> UnitCellObject;

    BravaisLatticeOfUnitCell(std::array<int, d> const& L, Real Scale = 1.):
      BravaisLatticeObject(L, Scale*ScaleBravais, AngleInGreekPis),
      UnitCellObject(Scale) {}
  };

  template<int d, char c, int i2, char c2>
  Real BravaisLatticeOfUnitCell<d, c, i2, c2>::ScaleBravais = 1.;

  template<int d, char c, int i2, char c2>
  Real BravaisLatticeOfUnitCell<d, c, i2, c2>::AngleInGreekPis = 0.;

}//namespace ffd::lattice
