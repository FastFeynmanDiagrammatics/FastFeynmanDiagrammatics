namespace ffd::principal_minors::unit_test{

  void matrix_6_by_6_step_2(){
	 using ffd::vector_range::Range;
     using ffd::set_theory::BinarySpread;

	 std::vector<Real> Mtx(36,0.0);
	 std::vector<Real> Rslt_minors(64,0.0);
	 std::vector<Real> Rslt_minors_step(8,0.0);
	 std::vector<Real> Rslt_minors_old(8,0.0);
	 std::vector<Real> Rslt_minors_2n(8,0.0);

	 for (int element : Range(Mtx.size())){
		 Mtx[element] = 2.0 * cos(element+1.0) -1.0;
 		 // Mtx[element] = 2.0 * ffd::random_distributions::Proba() -1.0;
	 }

	 Rslt_minors = PrincipalMinors(Mtx);
	 Rslt_minors_step = PrincipalMinors(Mtx,2);
	 Rslt_minors_old = PMOLD(Mtx,2);
	 Rslt_minors_2n = PrincipalMinors2n(Mtx,2);
	 for (int minor : Range(Rslt_minors_step.size())){
         // std::cerr << minor << " " << BinarySpread(minor, 2) << " " << Rslt_minors_step[minor] << " " << Rslt_minors_2n[minor] <<" " << Rslt_minors[BinarySpread(minor, 2)] << std::endl;
        assert(std::abs(Rslt_minors_step[minor] - Rslt_minors[BinarySpread(minor, 2)]) <1e-10);
	 }
  }

}//namespace
