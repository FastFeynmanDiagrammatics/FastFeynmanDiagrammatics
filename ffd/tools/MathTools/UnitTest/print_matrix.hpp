namespace ffd::math_tools::unit_test{

  void test_print_matrix(){
    std::vector<double> M(16,0.);

    for(uint j = 0; j < M.size(); ++j){
      M[j] = 2*ffd::random_distributions::Proba()-1;
    }

    print_matrix(M);

    print_matrix(M,4);

    print_matrix(M,4,4);

    print_matrix(M,2,8);

    for(uint j = 0; j< M.size(); ++j){
      M[j] *= 100;
    }

    print_matrix(M);

    print_matrix<10>(M);

    std::cerr<< 3.1652356 << std::endl;

    for(uint j = 0; j < 4; ++j){
      M[j+j*4] = 0;
    }
    M[0] = 1./0.;
    M[1] = -1./0.;
    M[2] = std::nan("1");
    print_matrix_signs(M);

  }

}//namespace
