namespace ffd::gauss_pfaffian::unit_test{

  template<int n = 4,
	   typename field_t = long double>
  struct matrix_antisym_test{
    std::array<field_t, n*(2*n-1)> el;

    matrix_antisym_test(std::array<field_t, n*(2*n-1)> A_): el(A_) {}
    
    field_t operator()(int j, int k) const{
      if( j > k ){
	return -(*this)(k, j);
      }else{
	return el[ lin(j, k, n) ];
      }
    }
  };


  
  template<typename field_t = long double>
  bool eight_by_eight(std::array<field_t, 28> A_array){
    bool IsOk = true;

    matrix_antisym_test<4, field_t> A(A_array);


    auto pf1 = [A](int x0, int x1){
		 return A(x0, x1);
	       };

    
    auto pf2 = [=](int x0, int x1, int x2, int x3){
		 return pf1(x0, x1) * A(x2, x3);
	       };

    auto PF2 = [=](int x0, int x1, int x2, int x3){
		 return
		   + pf2(x0, x1, x2, x3)
		   - pf2(x0, x2, x1, x3)
		   - pf2(x2, x1, x0, x3);
	       };
    
    
    auto pf3 = [=](int x0, int x1, int x2, int x3, int x4, int x5){
		 return PF2(x0, x1, x2, x3) * A(x4, x5);
		 };

    
    auto PF3 = [=](int x0, int x1, int x2, int x3, int x4, int x5){
		 return
		   + pf3(x0, x1, x2, x3, x4, x5)
		   - pf3(x0, x1, x2, x4, x3, x5)
		   - pf3(x0, x1, x4, x3, x2, x5)
		   - pf3(x0, x4, x2, x3, x1, x5)
		   - pf3(x4, x1, x2, x3, x0, x5);
		 };


    auto pf4 = [=](int x0, int x1, int x2, int x3, int x4, int x5, int x6, int x7){
		 return PF3(x0, x1, x2, x3, x4, x5) * A(x6, x7);
	       };


    auto PF4 = [=](int x0, int x1, int x2, int x3, int x4, int x5, int x6, int x7){
		 return
		   + pf4(x0, x1, x2, x3, x4, x5, x6, x7)
		   - pf4(x0, x1, x2, x3, x4, x6, x5, x7)
		   - pf4(x0, x1, x2, x3, x6, x5, x4, x7)
		   - pf4(x0, x1, x2, x6, x4, x5, x3, x7)
		   - pf4(x0, x1, x6, x3, x4, x5, x2, x7)
		   - pf4(x0, x6, x2, x3, x4, x5, x1, x7)
		   - pf4(x6, x1, x2, x3, x4, x5, x0, x7);
	       };

    
    auto pfaffian_exact = PF4(0, 1, 2, 3, 4, 5, 6, 7);
    auto pfaffian_gauss = Pfaffian(std::move(A_array));

    
    // std::cerr<<std::setprecision(20 )<<pfaffian_exact<<std::endl;
    // std::cerr<<pfaffian_gauss<<std::endl;

    
    IsOk = IsOk &&
      (
       std::abs(pfaffian_gauss - pfaffian_exact)
       < std::numeric_limits<field_t>::epsilon()
       *20*std::abs(pfaffian_exact)
       );

	    
    
    return IsOk;
  }

}//namespace
