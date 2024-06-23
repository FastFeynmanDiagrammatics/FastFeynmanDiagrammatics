

namespace ffd::lattice::unit_test{

  void creation_bravais_lattice_of_unit_cell(){
    using namespace std;

    BravaisLatticeOfUnitCell<2, '0', 2, '0'> B({2, 2});

    assert( abs( B.ScaleBravais - 1.) < numeric_limits<Real>::epsilon() );
    assert( abs( B.AngleInGreekPis ) < numeric_limits<Real>::epsilon() );


    HoneycombLattice H({4, 4});

    assert( abs( H.ScaleBravais - sqrt(3.)) < numeric_limits<Real>::epsilon() );
    assert( abs( H.AngleInGreekPis - 1./6) < numeric_limits<Real>::epsilon() );


    KagomeLattice K({5, 6});

    assert( abs( K.ScaleBravais - 2.) < numeric_limits<Real>::epsilon() );
    assert( abs( K.AngleInGreekPis) < numeric_limits<Real>::epsilon() );

  }

}//namespace
