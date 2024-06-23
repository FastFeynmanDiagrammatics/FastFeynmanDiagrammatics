namespace ffd::user_space{

  template<int d, typename... args>
  auto
  operator*(ffd::periodic_coordinate::PeriodicCoordinates<d, Real> const& tau_,
	    ffd::periodic_coordinate::PeriodicCoordinates<d, args...> const& x_){
    return ffd::merge(tau_, x_);
  }
  
}//namespace ffd::user_space
