namespace ffd::periodic_coordinate{
  
  template<typename T, int d>
  T PeriodicCoordinate<T, d>::operator()() const{
    if( !LowerUpperBound.has_value() ){
      return Variable;
    }
    auto const& [Lower, Upper] = LowerUpperBound.value();
    auto length = Upper - Lower;
    auto variable = Variable;
    while(variable < Lower){
      variable += length;
    }
    while(variable >= Upper){
      variable -= length;
    }
    return variable;
  }

}//namespace ffd::periodic_coordinate
