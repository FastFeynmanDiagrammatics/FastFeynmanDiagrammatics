

namespace ffd::principal_minors::unit_test{

  void matrix_4_by_4_step(){
	 using ffd::vector_range::Range;
	 using ffd::math_tools::print_matrix;
     using ffd::set_theory::BinarySpread;

	 std::vector<Real> Mtx(16,0.0);
	 std::vector<Real> Rslt_minors(16,0.0);
	 std::vector<Real> Rslt_minors_step(4,0.0);
 	 std::vector<Real> Rslt_minors_det(4,0.0);

	 for (int element : Range(Mtx.size())){
		 Mtx[element] = 2.0 * cos(element+1.0) -1.0;
 		 // Mtx[element] = 2.0 * ffd::random_distributions::Proba() -1.0;
	 }
  //   std::cerr << "4x4 " << std::endl;
	 Rslt_minors = PrincipalMinors(Mtx);
	 Rslt_minors_step = PrincipalMinors(Mtx,2);
	 Rslt_minors_det = DeterminantPrincipalMinors(Mtx,2);

	 for (int minor : Range(Rslt_minors_step.size())){
//          std::cerr << Rslt_minors_step[minor] << " " << Rslt_minors[BinarySpread(minor, 2)] << std::endl;
          assert(std::abs(Rslt_minors_step[minor] - Rslt_minors[BinarySpread(minor, 2)]) <1e-10);
		  assert(std::abs(Rslt_minors_step[minor] - Rslt_minors_det[minor]) <1e-10);
	 }
  }

}//namespace
