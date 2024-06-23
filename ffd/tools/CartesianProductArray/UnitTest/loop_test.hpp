

namespace ffd::cartesian_product_array::unit_test{

  void loop_test(){
    using ffd::vector_range::Range;

    CartesianProductArray<int, 3> r{{Range(2), Range(3), Range(5)}};

    int iter = 0;
    for(auto v: r){
      auto [x, y, z] = v;
      assert((
	      x == iter%2
	      ));
      assert((
	      y == (iter/2)%3
	      ));
      assert((
      	      z == (iter/2/3)%5
      	      ));
      ++iter;
    }


  }


}//namespace
