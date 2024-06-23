

namespace ffd::principal_minors::unit_test{

  void matrix_3_by_3_classic(){
	 using ffd::vector_range::Range;

	 std::vector<Real> Mtx(9,0.0);
	 std::vector<Real> Rslt(8,0.0);
	 std::vector<Real> Rslt_minors(8,0.0);
	 std::vector<Real> Rslt_dets(8,0.0);

	 for (int element : Range(Mtx.size())){
		 // Mtx[element] = 2.0 * ffd::random_distributions::Proba() -1.0;
     Mtx[element] = 2.0 * cos(element+1.0) -1.0;
	 }
   // Mtx = {4,7,3,6,2,9,6,1,8};
	 Rslt[0] = 1.0;
	 Rslt[1] = Mtx[0];
	 Rslt[2] = Mtx[4];
	 Rslt[4] = Mtx[8];
	 Rslt[3] = Mtx[0]*Mtx[4] - Mtx[1]*Mtx[3];
	 Rslt[5] = Mtx[0]*Mtx[8] - Mtx[2]*Mtx[6];
	 Rslt[6] = Mtx[4]*Mtx[8] - Mtx[5]*Mtx[7];

	 Rslt[7]  = Mtx[0]*Mtx[4]*Mtx[8];
	 Rslt[7] += Mtx[1]*Mtx[5]*Mtx[6];
	 Rslt[7] += Mtx[2]*Mtx[3]*Mtx[7];
	 Rslt[7] -= Mtx[2]*Mtx[4]*Mtx[6];
	 Rslt[7] -= Mtx[0]*Mtx[5]*Mtx[7];
	 Rslt[7] -= Mtx[1]*Mtx[3]*Mtx[8];

	 Rslt_minors = PrincipalMinorsBFS(Mtx);

	 for (int minor : Range(Rslt.size())){
		 // std::cerr << minor << " " << Rslt[minor] <<" " <<  Rslt_minors[minor] << std::endl;
	 	assert(std::abs(Rslt[minor] - Rslt_minors[minor]) <1e-10);
	 }

// exit(1); 

  }

}//namespace
