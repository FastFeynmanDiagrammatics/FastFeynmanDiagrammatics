namespace ffd::lattice{
  
  using TriangularLattice = BravaisLattice<2, 'T'>;

  void InitializeBravais(TriangularLattice& T_, Real lattice_constant){
    T_.BravaisVectors[0] = {lattice_constant, 0};
    Real x = lattice_constant*0.5l;
    Real y = lattice_constant*0.86602540378443864679l;
    T_.BravaisVectors[1] = {x, y};
  }

}//namespace ffd::lattice
