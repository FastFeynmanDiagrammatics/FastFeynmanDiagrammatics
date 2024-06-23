namespace ffd::gauss_pfaffian::unit_test{

  
  struct matrix_antisym6{
    std::array<long double, 15> el;

    matrix_antisym6(std::array<long double, 15> A_): el(A_) {}
    
    long double operator()(int j, int k) const{
      if( j > k ){
	return -(*this)(k, j);
      }else{
	return el[ lin(j, k, 3) ];
      }
    }
  };
  

  
  bool six_by_six(std::array<long double, 15> A_array){
    bool IsOk = true;


    matrix_antisym6 A(A_array);

    

    auto pf_ex = [A](int x0, int x1, int x2, int x3, int x4, int x5){
		   return ( + A(x0, x1) * A(x2, x3)
			    - A(x0, x2) * A(x1, x3)
			    + A(x0, x3) * A(x1, x2) ) * A(x4, x5);
		 };

    long double pfaffian_exact =
      +pf_ex(0, 1, 2, 3, 4, 5)
      -pf_ex(0, 1, 2, 4, 3, 5)
      -pf_ex(0, 1, 4, 3, 2, 5)
      -pf_ex(0, 4, 2, 3, 1, 5)
      -pf_ex(4, 1, 2, 3, 0, 5);
      
    auto pfaffian_gauss = Pfaffian(std::move(A_array));
    
    // std::cerr<<std::setprecision(20 )<<pfaffian_exact<<std::endl;
    // std::cerr<<pfaffian_gauss<<std::endl;

    
    IsOk = IsOk &&
      (std::abs(
		pfaffian_gauss - pfaffian_exact
		) < std::numeric_limits<long double>::epsilon()*20
       );

	    
    
    return IsOk;
  }

}//namespace



