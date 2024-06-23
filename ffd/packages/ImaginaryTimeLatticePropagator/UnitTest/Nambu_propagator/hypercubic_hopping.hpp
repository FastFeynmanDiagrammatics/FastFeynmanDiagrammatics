namespace ffd::user_space::imaginary_time_lattice_propagator::unit_test{


  template<int d, int j, typename T>

  auto

  hypercubic_hopping(T Origin,
		     std::array<std::array<Real, d>, 2> t){
    QuadraticAction<Real> S;

    
    for( int spin: {-1, 1} ){
      auto X = Origin;
      component<j+1>(X) = 1;
      auto hop = Bar(Psi_(spin))(X)*Psi_(spin)(Origin) * ( - t[(spin+1)/2][j] );
      S += hop + HermitianConjugate(hop);
    }

    
    if constexpr(j+1 < d){
	return S + hypercubic_hopping<d, j+1>(Origin, t);
      }else{
      return S;
    }
  }


}//namespace
