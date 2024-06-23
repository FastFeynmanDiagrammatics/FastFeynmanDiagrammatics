namespace ffd{

  template<int first, int... last, int d, typename... types>
  std::array<Real, d>
  RealSpaceFold(ffd::periodic_coordinate::PeriodicCoordinates<d, types...> const& PC_){
    using ffd::get;
    std::vector<std::array<Real, d>> VectorOfRet;
    ((VectorOfRet.push_back(ffd::RealSpace(get<first>(PC_)))), ... ,
     (VectorOfRet.push_back(ffd::RealSpace(get<last>(PC_)))));
    std::array<Real, d> ret;
    ret.fill(0.);
    for(auto vec: VectorOfRet){
      for(int j=0; j < d; ++j){
	ret[j] += vec[j];
      }
    }
    return ret;
  }

  
  template<int... last, int d, typename... types>
  std::array<Real, d>
  RealSpaceFold(ffd::periodic_coordinate::PeriodicCoordinates<d, types...> const& PC_,
		ffd::packages_math::StaticRangeSequence<0, 0, last...>){
    return RealSpaceFold<last...>(PC_);
  }

  
  template<int d, typename... types>
  std::array<Real, d>
  RealSpace(ffd::periodic_coordinate::PeriodicCoordinates<d, types...> const& PC_){
    ffd::packages_math::StaticRangeSequence<0, sizeof...(types)> S_;
    return RealSpaceFold(PC_, S_);
  }

  
}//namespace ffd
