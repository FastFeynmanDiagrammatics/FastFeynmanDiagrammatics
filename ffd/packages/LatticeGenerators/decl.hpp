namespace ffd::lattice_g{

  using int_d = int_fast16_t;
  
  template<int d, int n_atoms = 1>
  using coord_d = std::array<int_d, d+(n_atoms!=1)>;

}//namespace
