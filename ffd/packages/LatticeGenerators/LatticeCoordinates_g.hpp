namespace ffd::user_space{

  template<int d = 2,
	   int n_atoms = 1>

  auto
  LatticeCoordinates_g(){
    auto ret =
      []
      (auto... x){
	if constexpr(sizeof...(x) == d + (n_atoms!=1)){
	    return ffd::lattice_g::
	      coord_d<d, n_atoms>({static_cast<ffd::lattice_g::int_d>(x)...});
	  }else{
	  auto static_assert_error = [](){};
	  return static_assert_error;
	}
      };
    return ret;
  }

  
}//namespace
