

namespace ffd::eigenvalues_vectors::unit_test{

  void eigenvalues_vectors_call(){
    using ffd::user_space::I;
    

    std::vector<std::vector<Real>> S =
      {{
	{{1., 2.}},
	{{2., 2.}}
	}};



    auto [eigenvalues, eigenvectors] =
      EigenvaluesVectors<Real>(S);



    std::vector<std::vector<Complex>> H =
      {{
	{{1., 2. + I}},
	{{2. - I, 2.}}
	}};



    auto [eigenvalues_c, eigenvectors_c] =
      EigenvaluesVectors<Complex>(H);


  }


}//namespace
