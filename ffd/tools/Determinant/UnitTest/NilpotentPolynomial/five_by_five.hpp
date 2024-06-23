namespace ffd::determinant::unit_test{

#ifdef FFD_DETERMINANT_FIVE_UNIT_TEST_FLAG
  
  void five_by_five(){
    using nilpoly_t = ffd::nilpotent_polynomial::NilpotentPolynomial<Real>;
    using nilpolyl_t = ffd::nilpotent_polynomial::NilpotentPolynomial<long double>;
    long double const eps_real = std::numeric_limits<Real>::epsilon();
    int const matrix_linear_size = 5;
    std::vector<nilpoly_t> P(matrix_linear_size*matrix_linear_size);
    std::vector<nilpolyl_t> Pl(matrix_linear_size*matrix_linear_size);
    long double const omega = .96821821;
    long double const phi0 = .812312312432;
    int const order = 4;
    int const order_check = 1;


    for(int j=0; j<matrix_linear_size*matrix_linear_size; ++j){
      P[j] = nilpoly_t(order);
      Pl[j] = nilpolyl_t(order);
      for( BinaryInt S=0; S < (1<<order); ++S){
	long double pseudo_random = omega*(S+j*4+10) + phi0;
	pseudo_random = 2*(pseudo_random - std::round(pseudo_random));
	// std::cerr<<j<<" "<<S<<" "<<pseudo_random<<std::endl;
	P[j][S] = pseudo_random;
	Pl[j][S] = pseudo_random;
      }
    }
    
    
    auto det_poly = Determinant(P);

    
    auto vector_det = det5_test(Pl);


    for( BinaryInt S=0; S < (1<<order_check); ++S){
      std::cerr<<S<<" = "
      	       <<vector_det[S]<<" "
	       <<det_poly[S]<<" "
      	       // <<det_poly[S]-vector_det[S]<<" "
      	       <<std::endl;
      // assert((
      // 	      std::abs(vector_det[S] - det_poly[S]) <
      // 	      100*eps_real
      // 	      ));
    }
    
  }
  
#endif
  
}//namespace
