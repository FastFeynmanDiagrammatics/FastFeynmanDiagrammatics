

namespace ffd::principal_minors::unit_test{

  void matrix_4_by_4_Leading(){
	 using ffd::vector_range::Range;

   std::vector<Real> Mtx(6*6,0.0);
   std::vector<Real> Rslt_minors_A((1<<6),0.0);
   std::vector<Real> Rslt_minors_L(6,0.0);
   std::vector<Real> Rslt_minors_2L(3,0.0);
   std::vector<Real> Rslt_det_minors_L(6,0.0);
   std::vector<Real> Rslt_det_minors_2L(6,0.0);

	 for (int element : Range(Mtx.size())){
		 // Mtx[element] = ffd::random_distributions::Proba();
 		 Mtx[element] = 2.0 * ffd::random_distributions::Proba() -1.0;
	 }

   auto Mtx_ = Mtx;
   auto Mtx__ = Mtx;
   auto Mtx___ = Mtx;

	 Rslt_minors_L = LeadingPrincipalMinors(Mtx);
 	 Rslt_minors_2L = LeadingPrincipalMinors2n(Mtx_);
   Rslt_det_minors_L = DeterminantLeadingPrincipalMinors(Mtx__);
   Rslt_det_minors_2L = DeterminantLeadingPrincipalMinors(Mtx___,2);


	 for (int minor : Range(Rslt_minors_L.size())){
		 // std::cerr << minor << " " << Rslt_minors_DFS[minor] <<" " <<  Rslt_minors_DET[minor] << std::endl;
	 	assert(std::abs(Rslt_minors_L[minor] - Rslt_det_minors_L[minor]) <1e-10);
    assert(std::abs(Rslt_minors_2L[minor] - Rslt_det_minors_2L[minor]) <1e-10);
	 }

   // for (int minor : Range(Rslt_minors_DET.size())){
   //   std::cerr << minor << " "
   //             << Rslt_minors_DET[minor] << " "
   //             << Rslt_minors_DFS[minor] << " "
   //             << Rslt_minors_BFS[minor] << " "
   //             << std::endl;
   // }

  }



}//namespace
