namespace ffd::lattice_g{

  template<int d, int n_atoms = 1>

  auto
  Realspace_g(std::array<int, d> const& L,
	      std::array<std::array<Real, d>, d> const& bravais,
	      std::array<std::array<Real, d>, n_atoms> const& shift_atoms
	      = std::array<std::array<Real, d>, n_atoms>()){
    auto r1 = Realspace_unary_g<d, n_atoms>(L, bravais, shift_atoms);
    auto r2 = Realspace_diff_g<d, n_atoms>(L, bravais, shift_atoms);

    auto ret =
      [r1, r2]
      (auto&&... x){
	if constexpr(sizeof...(x) == 1){
	    return r1(std::forward<decltype(x)>(x)...);
	  }else{
	  return r2(std::forward<decltype(x)>(x)...);
	}
      };

    return ret;
  }
}//namespace
