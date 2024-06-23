

namespace ffd::packages_math::unit_test{

  void static_range_sequence(){
    static_assert(std::is_base_of<StaticRangeSequence<0, 0, 0>, StaticRangeSequence<0, 1>>::value);
    static_assert(std::is_base_of<StaticRangeSequence<3, 3, 3, 4, 5, 6>, StaticRangeSequence<3, 7>>::value);
    static_assert(std::is_base_of<StaticRangeSequence<0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9>,
		  StaticRangeSequence<0, 10>>::value);
    static_assert(std::is_base_of<StaticRangeSequence<-1, -1, -1, 0, 1>,
		  StaticRangeSequence<-1, 2>>::value);
	
  }

}//namespace
