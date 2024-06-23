namespace ffd::matrix_2x2::unit_test{

	using namespace ffd::principal_minors;

	void test_principal_minors(){

    using namespace ffd::user_space;

		std::vector<Matrix2x2> Matrix2(4);
    std::vector<Real> Matrix4(16);

    for (uint i=0; i<16; ++i){
      Matrix4[i] = Proba();
    }

    Matrix2[0].coeff[0] = Matrix4[0];
    Matrix2[0].coeff[1] = Matrix4[1];
    Matrix2[0].coeff[2] = Matrix4[4];
    Matrix2[0].coeff[3] = Matrix4[5];

    Matrix2[1].coeff[0] = Matrix4[2];
    Matrix2[1].coeff[1] = Matrix4[3];
    Matrix2[1].coeff[2] = Matrix4[6];
    Matrix2[1].coeff[3] = Matrix4[7];

    Matrix2[2].coeff[0] = Matrix4[8];
    Matrix2[2].coeff[1] = Matrix4[9];
    Matrix2[2].coeff[2] = Matrix4[12];
    Matrix2[2].coeff[3] = Matrix4[13];

    Matrix2[3].coeff[0] = Matrix4[10];
    Matrix2[3].coeff[1] = Matrix4[11];
    Matrix2[3].coeff[2] = Matrix4[14];
    Matrix2[3].coeff[3] = Matrix4[15];

    auto PM2 = PrincipalMinors(Matrix2);
    auto PM4 = PrincipalMinors(Matrix4);

    std::cout << "principal minors of Real matrix " << std::endl;
    for (uint i=0; i<(1<<4); ++i){
      std::cout << i << " " << PM4[i] << std::endl;
    }
    std::cout << "principal minors of nonsense matrix " << std::endl;
    for (uint i=0; i<(1<<2); ++i){
      std::cout << i << " " << PM2[i].det() << std::endl;
    }




	}

}//namespace
