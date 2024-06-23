

namespace ffd::principal_minors::unit_test{

  void test_bounded(){
	 using ffd::vector_range::Range;

   const BinaryInt order = 5;
   const BinaryInt bound = 3;

	 std::vector<Real> Mtx(order*order,0.0);
	 std::vector<Real> Rslt_minors_DFS((1<<order),0.0);
 	 // std::vector<Real> Rslt_minors_BFS((1<<6),0.0);
	 std::vector<Real> Rslt_minors_DET((1<<order),0.0);

	 for (int element : Range(Mtx.size())){
		 Mtx[element] = 2.0 * ffd::random_distributions::Proba() -1.0;
     // Mtx[element] = 2.0 * cos(element+1.0) -1.0;
	 }

   Rslt_minors_DET = DeterminantPrincipalMinors<bound>(Mtx);
 	 Rslt_minors_DFS = PrincipalMinorsBounded(Mtx,bound);


	 for (int minor : Range(Rslt_minors_DET.size())){
		 // std::cerr << minor
     //           << " " << __builtin_popcount(minor)
     //           << " " << Rslt_minors_DFS[minor]
     //           << " " << Rslt_minors_DET[minor]
     //           << std::endl;
	 	assert(std::abs(Rslt_minors_DFS[minor] - Rslt_minors_DET[minor]) <1e-10);
    if (__builtin_popcount(minor)>bound){
      assert(Rslt_minors_DFS[minor]==0);
    }
	 }

  }



}//namespace
