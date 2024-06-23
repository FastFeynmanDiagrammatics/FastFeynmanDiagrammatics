namespace ffd::inverse_matrix{

  template<typename Field>
  std::vector<Field>
  InverseOfMatrix(std::vector<Field> M_){
    using std::abs;
    std::vector<Field> x(size(M_), 0.);



    int const size = ffd::core_math::sqrt_int(std::size(M_));


    
    std::vector<Field> rhs(std::size(M_), 0.);
    for(int j = 0; j < size; ++j){
      rhs[j+size*j] = 1.;
    }



    for(int j = 0; j < size; ++j){
      Real pivot_abs = abs(M_[j+j*size]);
      int pivot_position = j;
      for(int k = j+1; k < size; ++k){
	if(  pivot_abs  <  abs( M_[j+k*size] )  ){
	  pivot_position = k;
	  pivot_abs = abs( M_[j+k*size] );
	}
      }

      

      if( pivot_position != j ){
	SwapRows(M_ , {j, pivot_position});
	SwapRows(rhs, {j, pivot_position});
      }



      Field pivot_value = M_[j+j*size];
      MultiplyRow1OfMatrix2By3(j, M_ , 1./pivot_value);
      MultiplyRow1OfMatrix2By3(j, rhs, 1./pivot_value);



      for(int k = j + 1; k < size; ++k){
	Field row_first_element = M_[j+k*size];

	
	Row1Plus2Row3Matrix4(k, -row_first_element, j, M_);
	Row1Plus2Row3Matrix4(k, -row_first_element, j, rhs);
      }
    }


    x = rhs;
    for(int k = 0; k < size; ++k){
      for(int j = size-2; j >= 0; --j){
	for(int m = j + 1; m < size; ++m){
	  x[size*j + k]  -=  M_[size*j + m] * x[size*m + k];
	}
      }
    }
    return x;
  }


}//namespace
