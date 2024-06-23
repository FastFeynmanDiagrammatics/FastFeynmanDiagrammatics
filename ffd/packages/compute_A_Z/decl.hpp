namespace ffd::user_space{

  template<int order,
	   typename field_t, 
	   typename nilpoly_t =
	   ffd::nilpotent_polynomial::
	   NilpotentPolynomial<field_t>>
  
  std::array<nilpoly_t, 2>
  
  compute_A_Z(std::vector<field_t> const& PM,
	      std::array<std::array<int, (1<<order)>, 2> const& signs){
    std::array<nilpoly_t, 2> ret;
    ret.fill(nilpoly_t(order));
    BinaryInt const shift_full =
      size(PM)-1 -((1<<order)-1);

    
    for(BinaryInt S=0; S<(1<<order); ++S){
      ret[0][S] = PM[S+shift_full]*Real(signs[0][S]);
      ret[1][S] = PM[S]*Real(signs[1][S]);
    }
    
    
    return ret;
  }
	
}//namespace
