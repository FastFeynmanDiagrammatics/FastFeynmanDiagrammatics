

namespace ffd::user_space::imaginary_time_lattice_propagator::unit_test{

  auto instantiation(){
    constexpr int d = 3;

    
    int L = 11;
    ffd::lattice::HypercubicLattice<d> Cubic{{L, L, L}};

    
    Real Beta = 2.28;
    ffd::imaginary_time::PeriodicImaginaryTime ImagTime(Beta);


    auto X = CreateCoordinates(ImagTime, Cubic);
    auto Y = X;


    Real t = -1.;
    QuadraticAction A = Psi_(1)(X)*Bar(Psi_(1)(Y))*t;
    component<1>(Y) += 2;
    A += Psi_(1)(X)*Bar(Psi_(1)(Y))*t;
    component<2>(Y) -= 1;
    A += Psi_(-1)(X)*Bar(Psi_(-1)(Y))*t;
    A += HermitianConjugate(A);


    
    
    ImaginaryTimeLatticePropagator G0(ImagTime,
				      Cubic,
				      A);

    Complex z = 1e-10;
    ImaginaryTimeLatticePropagator G02(ImagTime,
				       Cubic,
				       A,
				       z);
    
    
    return G0;

  }


}//namespace
