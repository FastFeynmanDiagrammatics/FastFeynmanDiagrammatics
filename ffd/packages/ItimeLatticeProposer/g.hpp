namespace ffd::itime_lattice_proposer{

  template<class time_proposer_f,
	   class lattice_proposer_f>
  auto
  g(time_proposer_f const& Time_Proposer,
    lattice_proposer_f const& Lattice_Proposer){
    auto ret =
      [Time_Proposer, Lattice_Proposer]
      (auto&&... x){
	if constexpr(sizeof...(x) == 1){
	    auto [time, space] =
	      std::get<0>(split_arguments(x...));
	    return std::make_pair(Time_Proposer(time),
				  Lattice_Proposer(space));
	  }else{
	  auto [X0, X1] = split_arguments(x...);
	  auto [t0, x0] = X0;
	  auto [t1, x1] = X1;
	  return Time_Proposer(t0, t1)*
	    Lattice_Proposer(x0, x1);
	}
      };
    
    return ret;
  }

}//namespace
