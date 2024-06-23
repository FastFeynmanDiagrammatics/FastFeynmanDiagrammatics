namespace ffd::truncated_polynomial::unit_test{

  void test_TP(){

    constexpr int order = 10;
    ffd::truncated_polynomial::P<order> P1, P2,P4,P5;
    for (uint ord=0; ord<=order; ++ord){
      P1[ord] = 2.*((double) rand() / (RAND_MAX))-1.;
      P2[ord] = 2.*((double) rand() / (RAND_MAX))-1.;
      P4[ord] = 2.*((double) rand() / (RAND_MAX))-1.;
      P5[ord] = 2.*((double) rand() / (RAND_MAX))-1.;
    }

    // std::cout << ToString(P1) << std::endl;
    // std::cout << ToString(P2) << std::endl;

    {
      auto P3 = 3.5*(P2/3.5);
      for (uint ord=0; ord<=order; ++ord){
        // std::cout << ord << " " << P1[ord] << " " << P2[ord] << " " << P3[ord] << std::endl;
        assert(std::abs(P3[ord]-P2[ord])<1e-10);
      }
    }

    {
      auto P3 = P2-P2;
      for (uint ord=0; ord<=order; ++ord){
        // std::cout << ord << " " << P2[ord] << " " << P3[ord] << std::endl;
        assert(std::abs(P3[ord]-0.)<1e-10);
      }
    }

    {
      auto P3 = -(P1/P2)*P2+P1;
      for (uint ord=0; ord<=order; ++ord){
        // std::cout << ord << " " << P1[ord] << " " << P2[ord] << " " << P3[ord] << std::endl;
        assert(std::abs(P3[ord]-0.)<1e-10);
      }
    }

    {
      auto P3 = (P1*P2)+(-P2*P1);
      for (uint ord=0; ord<=order; ++ord){
        // std::cout << ord << " " << P1[ord] << " " << P2[ord] << " " << P3[ord] << std::endl;
        assert(std::abs(P3[ord]-0.)<1e-10);
      }
    }

    {
      ffd::truncated_polynomial::P<order> P3 = 1.;
      P3+= -P3;
      for (uint ord=0; ord<=order; ++ord){
        // std::cout << ord << " " << P3[ord] << std::endl;
        assert(std::abs(P3[ord]-0.)<1e-10);
      }
    }

    {
      std::cout << "real testing starts here" << std::endl;
      std::array<ffd::truncated_polynomial::P<order>,4> mat;
      std::array<ffd::truncated_polynomial::P<order>,4> inv;
      std::array<ffd::truncated_polynomial::P<order>,4> res;
      std::array<ffd::truncated_polynomial::P<order>,4> res_;
      mat[0b00] = P1;
      mat[0b01] = P2;
      mat[0b10] = P4;
      mat[0b11] = P5;
      inv[0b00] =  mat[0b11];
      inv[0b01] = -mat[0b01];
      inv[0b10] = -mat[0b10];
      inv[0b11] =  mat[0b00];

      auto det = inv[0b00]*inv[0b11] - inv[0b01]*inv[0b10];
      inv[0b00] /= det;
      inv[0b01] /= det;
      inv[0b10] /= det;
      inv[0b11] /= det;

      res[0b00] = mat[0b00]*inv[0b00] + mat[0b01]*inv[0b10];
      res[0b01] = mat[0b00]*inv[0b01] + mat[0b01]*inv[0b11];
      res[0b10] = mat[0b10]*inv[0b00] + mat[0b11]*inv[0b10];
      res[0b11] = mat[0b10]*inv[0b01] + mat[0b11]*inv[0b11];

      res[0b00] -= 1.;
      res[0b11] -= 1.;

      for (uint ord=0; ord<=order; ++ord){
        std::cout << ord << " " <<res[0b00][ord] << " " << res[0b01][ord]<< " " << res[0b10][ord] << " " << res[0b11][ord]<< std::endl;
        assert(std::abs(res[0b00][ord]-0.)<1e-10);
        assert(std::abs(res[0b01][ord]-0.)<1e-10);
        assert(std::abs(res[0b10][ord]-0.)<1e-10);
        assert(std::abs(res[0b11][ord]-0.)<1e-10);
      }

      res_[0b00] = 1.;
      res_[0b11] = 1.;

      res_[0b00] -= mat[0b00]*inv[0b00] + mat[0b01]*inv[0b10];
      res_[0b01] -= mat[0b00]*inv[0b01] + mat[0b01]*inv[0b11];
      res_[0b10] -= mat[0b10]*inv[0b00] + mat[0b11]*inv[0b10];
      res_[0b11] -= mat[0b10]*inv[0b01] + mat[0b11]*inv[0b11];

      for (uint ord=0; ord<=order; ++ord){
        std::cout << ord << " " <<res_[0b00][ord] << " " << res_[0b01][ord]<< " " << res_[0b10][ord] << " " << res_[0b11][ord]<< std::endl;
        assert(std::abs(res_[0b00][ord]-0.)<1e-10);
        assert(std::abs(res_[0b01][ord]-0.)<1e-10);
        assert(std::abs(res_[0b10][ord]-0.)<1e-10);
        assert(std::abs(res_[0b11][ord]-0.)<1e-10);
      }
    }


    {
      // =,+,-,+=,*,*=,/=,unary-, =1.
      auto PS = shift(shift(P1, .5), -.5) - P1 ;
      for (uint ord=0; ord<=order; ++ord){
	//	std::cerr << ord << " " << PS[ord] << " " << PS[ord] << " " << PS[ord] << std::endl;
	assert(std::abs(PS[ord]-0.)<1e-10);
      }
    }

        {
      // =,+,-,+=,*,*=,/=,unary-, =1.
      auto PS = shift(P1, .5);
      //      std::cerr << P1(.5) << " " << PS(0.) << "\n";
      assert(std::abs(P1(.5)-PS(0.))<1e-10);
    }



  }



}
