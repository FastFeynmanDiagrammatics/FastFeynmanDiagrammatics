namespace ffd{

  template<typename T, int d>
  std::array<Real, d>
  RealSpace(ffd::periodic_coordinate::PeriodicCoordinate<T, d> const& x_){
    using std::size;
    
    std::array<Real, d> ret;
    ret.fill(0.);
    
    auto vec = x_.RealSpaceVectors;
    if(size(vec) == 1){
      for(int j=0; j < d; ++j){
	ret[j] = vec[0][j]*x_();
      }
    }else if(size(vec) > 1){
      for(int j=0; j < d; ++j){
	const int index = x_();
	ret[j] = vec[index][j];
      }
    }
    
    return ret;
  }


}//namespace ffd
