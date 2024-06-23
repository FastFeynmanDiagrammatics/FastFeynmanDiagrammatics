

namespace ffd::eigenvalues_vectors::unit_test{

  void two_by_two_hermitian_matrix(){
    using ffd::user_space::I;

    std::array<std::array<Complex, 2>, 2> H;
    H = {{
	  {{1.  ,  1.+I}},
	  {{1.-I,  1.  }}
	  }};
    assert((two_by_two_hermitian_matrix_test_arg(H)));
    
  }


}//namespace
