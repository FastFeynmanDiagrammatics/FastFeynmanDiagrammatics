namespace ffd::lattice{

  
  template<int d, char c, typename... args>
  auto
  CreateCoordinatesIncremental(ffd::lattice::BravaisLattice<d, c> const& L_,
			       ffd::periodic_coordinate::PeriodicCoordinates<d, args...>& X_old){
    if constexpr (d == 0){
	return X_old;
      }else{
      using ffd::get, ffd::merge;
      ffd::periodic_coordinate::PeriodicCoordinates<d, int> X;
      get<0>(X).LowerUpperBound = L_.LowerUpperBounds[sizeof...(args)];
      get<0>(X).RealSpaceVectors.resize(1);
      get<0>(X).RealSpaceVectors[0] = L_.BravaisVectors[sizeof...(args)];
      get<0>(X).Variable = 0;
      auto Merged = merge(X_old, X);
      if constexpr( sizeof...(args) + 1 >= d){
	  return Merged;
	}else{
	return CreateCoordinatesIncremental(L_, Merged);
      }
    }
  }

  

  template<int d, char c>
  auto CreateCoordinates(ffd::lattice::BravaisLattice<d, c> const& L_){
    ffd::periodic_coordinate::PeriodicCoordinates<d> temp;
    return CreateCoordinatesIncremental(L_, temp);
  }

  
}//namespace 
