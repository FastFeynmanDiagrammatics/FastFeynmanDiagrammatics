

namespace ffd::user_space::imaginary_time_lattice_propagator::unit_test{

  void action_conversion(){

    int linear_system_size = 10;
    ffd::lattice::HoneycombLattice H{{linear_system_size,
				      linear_system_size}};

    
    Real Beta = 2.;
    ffd::imaginary_time::PeriodicImaginaryTime ImagTime(Beta);

    
    auto X = CreateCoordinates(ImagTime, H);

    
    constexpr int space_dimensions = 2;
    constexpr bool has_unit_cell = true;
    using coordinates_type = typename
      ffd::imaginary_time_lattice::ImaginaryTimeLatticeCoordinates<space_dimensions,
								   has_unit_cell>::type;

    
    static_assert(
		  std::is_same<decltype(X),
		  coordinates_type>::value
		  );

    
  
    
    
    Real t = 1.;

    auto term = Psi_(1)(X)*Bar(Psi_(1))(X)*t;

    
    static_assert(
		  std::is_same<decltype(term),
		  QuadraticActionTerm<Real>>::value);

    
    auto Y = std::any_cast<coordinates_type>(term[0].Position);

    component<1>(Y) += 1;


    QuadraticAction A = term;
    static_assert( std::is_same<decltype(A),
		   QuadraticAction<Real>>::value );

    
    A += Psi_(-1)(X)*Bar(Psi_(-1)(Y))*t;


    A += HermitianConjugate(A);

    
  }
  
}//namespace
