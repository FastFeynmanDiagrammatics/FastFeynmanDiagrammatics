

namespace ffd::principal_minors::unit_test{

  void matrix_4_by_4_step_0_diag(){
	 using ffd::vector_range::Range;
     using ffd::set_theory::BinarySpread;

	 std::vector<Real> Mtx(16,0.0);
	 std::vector<Real> Rslt_minors(16,0.0);
	 std::vector<Real> Rslt_minors_step(4,0.0);

    // for (uint cnt=0; cnt<(1<<20); ++cnt){
	// if (cnt%10000==0) std::cerr << cnt << std::endl;
	 for (int element : Range(Mtx.size())){
		 Mtx[element] = (2.0 * cos(element+1.0) -1.0)/1.e10;
 		 // Mtx[element] = (2.0 * ffd::random_distributions::Proba() -1.0)/1.e10;
	 }
	 for (int el =0; el<4; ++el){
		 Mtx[el*4+el] = 0;
	 }

	 Rslt_minors = PrincipalMinors(Mtx);
	 Rslt_minors_step = PrincipalMinors(Mtx,2);
	 for (int minor : Range(Rslt_minors_step.size())){
          // std::cerr << minor << " " << BinarySpread(minor, 2) << " " << Rslt_minors_step[minor] << " " << Rslt_minors[BinarySpread(minor, 2)] << std::endl;
          assert(std::abs(Rslt_minors_step[minor] - Rslt_minors[BinarySpread(minor, 2)]) <1e-10);
	 }
 	// }
	// std::cerr << "done" << std::endl;
  }

}//namespace
