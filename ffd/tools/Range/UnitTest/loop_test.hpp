

namespace ffd::vector_range::unit_test{

  void loop_test(){


    int iter = 0;
    for(auto j: Range(-1, 2)){
      assert((
	      j == iter-1
	      ));
      ++iter;
    }
    
    assert((
	    iter == 3
	    ));

  }
  

}//namespace
