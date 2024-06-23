namespace ffd::determinant::unit_test{

  void two_by_two(){
    using nilpoly_t = ffd::nilpotent_polynomial::NilpotentPolynomial<Real>;
    Real const eps_real = std::numeric_limits<Real>::epsilon();
    std::array<nilpoly_t, 4> P;
    Real const omega = .723821821;
    Real const phi0 = .912312312432;
    int const order = 4;


    
    for(int j=0; j<4; ++j){
      P[j] = nilpoly_t(order);
      for( BinaryInt S=0; S < (1<<order); ++S){
	Real pseudo_random = omega*(S+j+10) + phi0;
	pseudo_random = 2*(pseudo_random - long(pseudo_random)) - 1.;
	P[j][S] = pseudo_random;
      }
    }
    

    auto vector_det = det2_test(P);


    auto det_poly = Determinant(P);


    for( BinaryInt S=0; S < (1<<order); ++S){
      // std::cerr<<S<<" "<<vector_det[S]<<" "<<det_poly[S]-vector_det[S]<<std::endl;
      assert((
	      std::abs(vector_det[S] - det_poly[S]) <
	      100*eps_real
	      ));
    }
    
  }

}//namespace
