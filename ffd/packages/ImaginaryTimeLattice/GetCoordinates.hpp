namespace ffd::imaginary_time_lattice{

  template<int j,
	   int d_eff,
	   typename coordinates_t>
  
  auto
  space_components_rec(coordinates_t X, std::array<int, d_eff> components){
    components[j] = ffd::get<j+1>(X)();
    if constexpr( j == d_eff - 1){
	return components;
      }else{
      return space_components_rec<j+1, d_eff>(X, components);
    }
  }
  


  template<typename lattice_t,
	   typename coordinates_t>
    
  auto
  GetSpaceCoordinates(coordinates_t X){
    constexpr int d_eff = lattice_t::dimension + (lattice_t::number_atoms_unit_cell>1);
    std::array<int, d_eff> components;

    
    return space_components_rec<0, d_eff>(X, components);
  }


  

  template<typename imaginary_time_t,
	   typename lattice_t,
	   typename coordinates_t>
  auto
  GetCoordinates(coordinates_t X){
    Real tau = ffd::get<0>(X)();
    auto space_x = GetSpaceCoordinates<lattice_t>(X);
    return std::make_pair(tau, space_x);
  }
  

}//namespace
