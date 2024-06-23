namespace ffd::lattice{

  using HoneycombLattice = BravaisLatticeOfUnitCell<2, 'T', 2, 'S'>;

  template<>
  Real HoneycombLattice::ScaleBravais = 1.7320508075688772935274463415058723l;

  template<>
  Real HoneycombLattice::AngleInGreekPis = .16666666666666666666666666666l;

}//namespace ffd::lattice
