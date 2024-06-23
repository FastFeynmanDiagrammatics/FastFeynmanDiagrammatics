namespace ffd::lattice{

  template<int d, char bravais_name, int n_atoms, char cell_name>
  auto CreateCoordinates(ffd::lattice::BravaisLatticeOfUnitCell
			 <d, bravais_name, n_atoms, cell_name> const& B_){
    using ffd::merge;
    auto Bravais_coord = CreateCoordinates(B_.BravaisLatticeObject);
    auto Cell_coord = CreateCoordinates(B_.UnitCellObject);
    return merge(Bravais_coord, Cell_coord);
  }

}//namespace ffd::lattice

