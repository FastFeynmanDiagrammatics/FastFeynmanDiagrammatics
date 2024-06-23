namespace ffd::user_space{

  template<typename value_t,
	   std::size_t size_v>
  
  struct SafeArray{
    std::array<value_t, size_v> Array;

    SafeArray() {}

    
    SafeArray(std::array<value_t, size_v> Array_):
      Array(Array_) {}


    value_t& operator[](std::size_t j){
      assert(( j >= 0 ));
      assert(( j < size_v));
      return Array[j];
    }
    
  };

  

}//namespace

