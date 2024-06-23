namespace ffd::principal_minors::unit_test{

  void matrix_6_by_6_step(){
	 using ffd::vector_range::Range;
     using ffd::set_theory::BinarySpread;

	 std::vector<Real> Mtx(36,0.0);
	 std::vector<Real> Rslt_minors(64,0.0);
	 std::vector<Real> Rslt_minors_step(4,0.0);

	 for (int element : Range(Mtx.size())){
		 // Mtx[element] = 2.0 * ffd::random_distributions::Proba() -1.0;
 		 Mtx[element] = 2.0 * cos(element+1.0) -1.0;
	 }

	 Rslt_minors = PrincipalMinors(Mtx);
	 Rslt_minors_step = PrincipalMinors(Mtx,3);
	 for (int minor : Range(Rslt_minors_step.size())){
//          std::cerr << minor << " " << BinarySpread(minor, 3) << " " << Rslt_minors_step[minor] << " " << Rslt_minors[BinarySpread(minor, 3)] << std::endl;
          assert(std::abs(Rslt_minors_step[minor] - Rslt_minors[BinarySpread(minor, 3)]) <1e-10);
	 }
  }

}//namespace
