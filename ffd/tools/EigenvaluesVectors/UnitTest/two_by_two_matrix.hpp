

namespace ffd::eigenvalues_vectors::unit_test{


  void two_by_two_matrix(){

    
    array<array<Real, 2>, 2> S;
    S = {{
	  {{1., 0.}},
	  {{0., 1.}}
	  }};
    assert((get<0>(two_by_two_matrix_test_arg(S))));


    
    S = {{
	  {{2., .2}},
	  {{.2, 1.}}
	  }};
    assert((get<0>(two_by_two_matrix_test_arg(S))));


    

    S = {{
	  {{1., .2}},
	  {{.2, 1.}}
	  }};
    assert((get<0>(two_by_two_matrix_test_arg(S))));


    
    S = {{
	  {{1., .2}},
	  {{.2, 2.}}
	  }};
    assert((get<0>(two_by_two_matrix_test_arg(S))));

    

    S = {{
	  {{1., 2}},
	  {{2, 2.}}
	  }};
    assert((get<0>(two_by_two_matrix_test_arg(S))));


    
    S = {{
	  {{1., -2}},
	  {{-2, 2.}}
	  }};
    assert((get<0>(two_by_two_matrix_test_arg(S))));



    S = {{
	  {{4., -2}},
	  {{-2, 2.}}
	  }};
    assert((get<0>(two_by_two_matrix_test_arg(S))));


    S = {{
	  {{4., 0}},
	  {{0, 2.}}
	  }};
    assert((get<0>(two_by_two_matrix_test_arg(S))));


    
    S = {{
	  {{0., 1}},
	  {{1, 0.}}
	  }};
    assert((get<0>(two_by_two_matrix_test_arg(S))));



    S = {{
	  {{0., 1}},
	  {{1, 1.}}
	  }};
    assert((get<0>(two_by_two_matrix_test_arg(S))));


    S = {{
	  {{1., 1}},
	  {{1, 1.}}
	  }};
    assert((get<0>(two_by_two_matrix_test_arg(S))));


  }
  

}//namespace
