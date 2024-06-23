namespace ffd::lattice{

  template<int space_dimensions, char bravais_lattice_name = '0'>
  class BravaisLattice{
  public:
    static constexpr int dimension = space_dimensions;

    static constexpr int number_atoms_unit_cell = 1;
    
    std::array<std::array<Real, space_dimensions>, space_dimensions> BravaisVectors;

    std::array<std::array<int, 2>, space_dimensions> LowerUpperBounds;
    
    BravaisLattice(std::array<int, space_dimensions> const& Linear_system_sizes,
		   Real scale_factor_of_the_lattice = 1.,
		   Real angle_of_XY_rotation_in_greek_pis = 0.);
  };
  
  
}//namespace ffd::lattice


