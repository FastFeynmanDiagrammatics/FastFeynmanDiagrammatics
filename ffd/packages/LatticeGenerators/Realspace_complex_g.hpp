namespace ffd::lattice_g{

  template<int d, int n_atoms = 1>

  auto
  Realspace_complex_g(std::array<int, d> const& L,
		      std::array<std::array<Real, d>, d> const& bravais,
		      std::array<std::array<Real, d>, n_atoms> const& shift_atoms
		      = std::array<std::array<Real, d>, n_atoms>()){
    std::array<Real, d> Pi2_L, L_Pi2;
    for(uint j=0; j<d; ++j){
      Pi2_L[j] = 2*ffd::core_math::Pi/L[j];
      L_Pi2[j] = 1./Pi2_L[j];
    }
    auto Realspace =
      [bravais, shift_atoms, Pi2_L, L_Pi2]
      (coord_d<d, n_atoms> const& x){
	std::array<Complex, d> r;
	r.fill(0.);
	for(uint j=0; j<d; ++j){
	  Complex const z = std::exp(Complex(0., Pi2_L[j]*x[j]))*L_Pi2[j];
	  for(uint l=0; l<d; ++l){
	    r[l] += bravais[j][l]*z;
	  }
	}
	if constexpr(n_atoms != 1){
	    for(uint j=0; j<d; ++j){
	      r[j] += Complex(0., shift_atoms[x[d]][j]);
	    }
	  }
	return r;
      };
    return Realspace;
  }


}//namespace
