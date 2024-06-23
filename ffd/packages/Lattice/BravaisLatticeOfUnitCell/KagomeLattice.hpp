

namespace ffd::lattice{

  //a Kagome is a triangular bravais lattice of equilater triangles
  using KagomeLattice = BravaisLatticeOfUnitCell<2, 'T', 3, 'E'>;

  template<>
  Real KagomeLattice::ScaleBravais = 2.;

}//namespace
