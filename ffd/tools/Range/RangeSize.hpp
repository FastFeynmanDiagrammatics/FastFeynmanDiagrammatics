namespace ffd::vector_range{

  template<typename sizeable_t>
  auto
  RangeSize(sizeable_t x){
    return Range( size(x) );
  }

  
}//namespace
