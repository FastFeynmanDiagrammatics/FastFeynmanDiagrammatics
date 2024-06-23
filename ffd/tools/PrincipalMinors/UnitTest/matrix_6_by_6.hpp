

namespace ffd::principal_minors::unit_test{

  void matrix_6_by_6_Real(){
	 using ffd::vector_range::Range;

	 std::vector<Real> Mtx(6*6,0.0);
	 std::vector<Real> Rslt_minors_DFS((1<<6),0.0);
 	 std::vector<Real> Rslt_minors_BFS((1<<6),0.0);
	 std::vector<Real> Rslt_minors_DET((1<<6),0.0);

 	 std::vector<Real> Rslt_minors_2n_DET((1<<3),0.0);
	 std::vector<Real> Rslt_minors_2n_BFS((1<<3),0.0);

	 for (int element : Range(Mtx.size())){
	    Mtx[element] = 2.0 * ffd::random_distributions::Proba() -1.0;
	 }

   Rslt_minors_DET = DeterminantPrincipalMinors(Mtx);
 	 Rslt_minors_DFS = PrincipalMinors(Mtx);
	 Rslt_minors_BFS = PrincipalMinorsBFS(Mtx);

   Rslt_minors_2n_DET = PrincipalMinors2n(Mtx);
  // Rslt_minors_2n_DET = DeterminantPrincipalMinors(Mtx,2);
   Rslt_minors_2n_BFS = PrincipalMinorsBFS2n(Mtx);

   for (int minor : Range(Rslt_minors_2n_DET.size())){
    // std::cerr << minor << " " << Rslt_minors_2n_BFS[minor] <<" " <<  Rslt_minors_2n_DET[minor] << std::endl;
   assert(std::abs(Rslt_minors_BFS[minor] - Rslt_minors_DET[minor]) <1e-10);
  }

	 for (int minor : Range(Rslt_minors_DET.size())){
		 // std::cerr << minor << " " << Rslt_minors_DFS[minor] <<" " <<  Rslt_minors_DET[minor] << std::endl;
	 	assert(std::abs(Rslt_minors_DFS[minor] - Rslt_minors_DET[minor]) <1e-10);
	 }


	 for (int minor : Range(Rslt_minors_DET.size())){
		 // std::cerr << minor << " " << Rslt_minors_BFS[minor] <<" " <<  Rslt_minors_DET[minor] << std::endl;
	 	assert(std::abs(Rslt_minors_BFS[minor] - Rslt_minors_DET[minor]) <1e-10);
	 }

  }



}//namespace
