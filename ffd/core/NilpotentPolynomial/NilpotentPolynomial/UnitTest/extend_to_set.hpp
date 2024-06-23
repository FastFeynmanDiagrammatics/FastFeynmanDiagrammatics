namespace ffd::nilpotent_polynomial::unit_test{

  void extend_to_set(){
    
    using nilpoly_t = NilpotentPolynomial<Real>;
    Real const eps_real = std::numeric_limits<Real>::epsilon();
    using std::abs;
    

    {
    nilpoly_t P(1);
    P[0] = 1.;
    P[1] = 2.;
    auto Q = P.extend_to_set(2, 2);
    

    assert(( abs(P[0]-Q[0])<
	     eps_real ));

    assert(( abs(Q[1])<
	     eps_real ));


    assert(( abs(P[1]-Q[2])<
	     eps_real ));

    
    assert(( abs(Q[3])<
	     eps_real ));
    }


    {
      nilpoly_t P(2);
      for(int j=0; j< (1<<2); ++j ){
	P[j] = j+1;
      }


      auto Q = P.extend_to_set(2 + 8, 9);

      
      assert(( abs(Q[0]-P[0])<
	       eps_real ));


      assert(( abs(Q[1]) < eps_real ));


      assert(( abs(Q[2]-P[1]) < eps_real ));


      assert(( abs(Q[3]) < eps_real ));


      assert(( abs(Q[4]) < eps_real ));


      assert(( abs(Q[5]) < eps_real ));


      assert(( abs(Q[6]) < eps_real ));


      assert(( abs(Q[7]) < eps_real ));


      assert(( abs(Q[8]-P[2]) < eps_real ));


      assert(( abs(Q[9]) < eps_real ));


      assert(( abs(Q[10]-P[3]) < eps_real ));


      assert(( abs(Q[11]) < eps_real ));


      assert(( abs(Q[12]) < eps_real ));


      assert(( abs(Q[13]) < eps_real ));


      assert(( abs(Q[14]) < eps_real ));

      
      assert(( abs(Q[15]) < eps_real ));

      
      assert(( abs(Q[16]) < eps_real ));

      
      assert(( abs(Q[17]) < eps_real ));
      
    }




    
    {

      
      nilpoly_t P(3);
      for( int j=0; j < 8; ++j){
	P[j] = j+1;
      }


      auto Q = P.extend_to_set(1+4+8, 6);


      assert(( abs(P[0]-Q[0]) < eps_real ));


      assert(( abs(P[1]-Q[1]) < eps_real ));


      assert(( abs(Q[2]) < eps_real ));


      assert(( abs(Q[3]) < eps_real ));


      assert(( abs(Q[4]-P[2]) < eps_real ));


      assert(( abs(Q[5]-P[3]) < eps_real ));


      assert(( abs(Q[6]) < eps_real ));


      assert(( abs(Q[7]) < eps_real ));


      assert(( abs(Q[8]-P[4]) < eps_real ));


      assert(( abs(Q[9]-P[5]) < eps_real ));


      assert(( abs(Q[10]) < eps_real ));


      assert(( abs(Q[11]) < eps_real ));


      assert(( abs(Q[12]-P[6]) < eps_real ));


      assert(( abs(Q[13]-P[7]) < eps_real ));


      assert(( abs(Q[14]) < eps_real ));


      assert(( abs(Q[15]) < eps_real ));


      assert(( abs(Q[15]) < eps_real ));
            

    }
    
  }

}//namespace
