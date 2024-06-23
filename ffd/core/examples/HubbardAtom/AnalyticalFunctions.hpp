namespace ffd::user_space::core_examples{

  auto EnergyOneSite = [](Real mu0, Real alpha, Complex U, std::pair<int, int> n)->Complex{
			 auto [n_down, n_up] = n;
			 return U*(n_up - alpha)*(n_down - alpha) - mu0*(n_up+n_down);
		       };

  
  std::function<Complex(Complex)> BuildPartitionFunction(Real InverseTemperature, Real ChemicalPotential){
    return [beta = InverseTemperature,
	    mu = ChemicalPotential,
	    alpha = 1./(1+exp(-InverseTemperature*ChemicalPotential))]
      (Complex U)->Complex{
	     Complex ret = 0;
	     for(auto n_up: {0, 1}){
	       for(auto n_down: {0, 1}){
		 ret += exp( - EnergyOneSite(mu, alpha, U, {n_down, n_up}) * beta );
	       }
	     }
	     return ret;;
	   };
  }

}  
  std::function<Complex(Complex)> BuildDensity(Real InverseTemperature, Real ChemicalPotential){
    auto Z_func = BuildPartitionFunction(InverseTemperature, ChemicalPotential);
    auto A_func = [beta = InverseTemperature,
	    mu = ChemicalPotential,
	    alpha = 1./(1+exp(-InverseTemperature*ChemicalPotential))]
      (Complex U)->Complex{
	     Complex ret = 0;
	     for(auto n_up: {0, 1}){
	       for(auto n_down: {0, 1}){
		 ret += exp( - EnergyOneSite(mu, alpha, U, {n_down, n_up}) * beta )*.5*(n_up+n_down+0.);
	       }
	     }
	     return ret;
	   };
    return [A_func, Z_func](Complex U)->Complex{return A_func(U)/Z_func(U);};
  }


  std::function<Complex(Complex)> BuildDoubleOccupancy(Real InverseTemperature, Real ChemicalPotential){
    auto Z_func = BuildPartitionFunction(InverseTemperature, ChemicalPotential);
    auto A_func = [beta = InverseTemperature,
		   mu = ChemicalPotential,
		   alpha = 1./(1+exp(-InverseTemperature*ChemicalPotential))]
      (Complex U)->Complex{
		    Complex ret = 0;
		    for(auto n_up: {0, 1}){
		      for(auto n_down: {0, 1}){
			ret += exp( - EnergyOneSite(mu, alpha, U, {n_down, n_up}) * beta )*(n_up-alpha)*(n_down-alpha);
		      }
		    }
		    return ret;
		  };
    
    return [A_func, Z_func](Complex U)->Complex{return A_func(U)/Z_func(U);};
  }
  
}//namespace
