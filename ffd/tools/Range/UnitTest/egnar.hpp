namespace ffd::vector_range::unit_test{


  void egnar(){
    
    int counter = 0;
    for( int j: egnaR(3) ){
      // std::cerr<<j<<std::endl;
      if(counter == 0){
	assert(( j == 2 ));
      }else if(counter == 1){
	assert(( j == 1 ));
      }else{
	assert(( j==0 ));
      }
      ++counter;
    }
    assert(( counter == 3 ));
    
  }

}//namespace
