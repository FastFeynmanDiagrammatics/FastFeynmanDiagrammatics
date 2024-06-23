namespace ffd::user_space{

  template<int d = 2,
	   int n_atoms = 1>
  
  auto
  Spacetime_g(){
    auto ret =
      []
      (Real x0, auto... x){
	static_assert(sizeof...(x) == d+(n_atoms!=1),
		      "check the number of parameters!!!");
	return std::make_pair(x0,
			      ffd::lattice_g::
			      coord_d<d, n_atoms>({static_cast<ffd::lattice_g::int_d>(x)...}));
      };
    return ret;
  }

}//namespace
