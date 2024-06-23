namespace ffd::nilpotent_polynomial::unit_test{

  void shift(){
    using std::abs;
    using nilpoly_t = NilpotentPolynomial<Real>;
    Real const eps_real = std::numeric_limits<Real>::epsilon();


    {
      nilpoly_t P(0);
      P[0] = 1.;
      auto Q = P.shift(1, 3);

      
      assert(( abs(Q[0]) < eps_real ));

      
      assert(( abs(Q[1]-P[0]) < eps_real ));


      assert(( abs(Q[2]) < eps_real ));


      assert(( abs(Q[3]) < eps_real ));

      
      assert(( abs(Q[4]) < eps_real ));
    }


    {
      nilpoly_t P(1);
      P[0] = 1.;
      P[1] = 2.;


      auto Q = P.shift(1+2, 3);


      assert(( abs(Q[0]) < eps_real ));


      assert(( abs(Q[1]) < eps_real ));


      assert(( abs(Q[2]) < eps_real ));


      assert(( abs(Q[3]-P[0]) < eps_real ));


      assert(( abs(Q[4]) < eps_real ));


      assert(( abs(Q[5]) < eps_real ));


      assert(( abs(Q[6]) < eps_real ));


      assert(( abs(Q[7]) < eps_real ));

    }

    
    {
      nilpoly_t P(3);
      for(int j=0; j< (1<<3); ++j){
	P[j] = j+1;
      }


      auto Q = P.shift(2, 3);


      assert(( abs(Q[0]) < eps_real ));


      assert(( abs(Q[1]) < eps_real ));


      assert(( abs(Q[2]-P[0]) < eps_real ));


      assert(( abs(Q[3] - P[1]) < eps_real ));

      
      assert(( abs(Q[4]) < eps_real ));


      assert(( abs(Q[5]) < eps_real ));


      assert(( abs(Q[6]-P[4]) < eps_real ));


      assert(( abs(Q[7] - P[5]) < eps_real )); 

    }


    {
      nilpoly_t P(4);
      for( int j=0; j < (1<<4); ++j){
	P[j] = j+1;
      }


      auto Q = P.shift(2+4, 4);


      assert(( abs(Q[0]) < eps_real ));


      assert(( abs(Q[1]) < eps_real ));


      assert(( abs(Q[2]) < eps_real ));

      
      assert(( abs(Q[3]) < eps_real ));


      assert(( abs(Q[4]) < eps_real ));


      assert(( abs(Q[5]) < eps_real ));


      assert(( abs(Q[6]-P[0]) < eps_real ));


      assert(( abs(Q[7] - P[1]) < eps_real ));


      assert(( abs(Q[8]) < eps_real ));


      assert(( abs(Q[9]) < eps_real ));


      assert(( abs(Q[10]) < eps_real ));


      assert(( abs(Q[11]) < eps_real ));


      assert(( abs(Q[12]) < eps_real ));

      
      assert(( abs(Q[13]) < eps_real ));


      assert(( abs(Q[14] - P[8]) < eps_real ));


      assert(( abs(Q[15] - P[9]) < eps_real ));

    }
  }

}//namespace
