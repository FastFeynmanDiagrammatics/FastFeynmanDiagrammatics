#pragma once

namespace ffd::cartesian_product::unit_test{

  void vector_range(){

    std::vector<int> x(1);
    x[0] = 1;
    for( auto r: ffd::user_space::VectorRange<1>(x) ){
      assert(( std::size(r) == 1));
      assert(( r[0] == 0 ));
    }

    

    std::vector<int> y(2);
    y[0] = 3; y[1] = 2;

    int counter = 0;
    for( auto r: ffd::user_space::VectorRange<2>(y) ){
      if(counter==0){
	assert(( r[0] == 0 ));
	assert(( r[1] == 0 ));
      }else if(counter==1){
	assert(( r[0] == 0 ));
	assert(( r[1] == 1 ));
      }else if(counter==2){
	assert(( r[0] == 1 ));
	assert(( r[1] == 0 ));
      }else if(counter==3){
	assert(( r[0] == 1 ));
	assert(( r[1] == 1 ));
      }else if(counter==4){
	assert(( r[0] == 2 ));
	assert(( r[1] == 0 ));
      }else if(counter==5){
	assert(( r[0] == 2 ));
	assert(( r[1] == 1 ));
      }
      ++counter;
    }
    assert(counter == 6);
   
    
  }


}//namespace
