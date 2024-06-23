

namespace ffd::principal_minors_old::unit_test{

  void matrix_3_by_3_Poly(){

  	using namespace ffd::nilpotent_polynomial;

  	using ffd::vector_range::Range;



   // NilpotentPolynomial<Real> P1(1);
   // P1[0] = 3.0;
   // P1[1] = 4.0;
   // NilpotentPolynomial<Real> P2 = 1.0/P1;
   //
   // std::cerr << P2[0] << " " << P2[1] << std::endl;
   // std::cerr << "popcnt = " << P2.size() << std::endl;
   // exit(1);

    int lin_size = 3;
    int cardinality = (1<<lin_size);

  	std::vector<NilpotentPolynomial<Real>> Mtx(lin_size*lin_size,NilpotentPolynomial<Real>(lin_size));
  	std::vector<NilpotentPolynomial<Real>> Rslt(cardinality,NilpotentPolynomial<Real>(lin_size));
  	std::vector<NilpotentPolynomial<Real>> Rslt_minors(cardinality,NilpotentPolynomial<Real>(lin_size));

  	for (int element : Range(Mtx.size())){
      for (int set : Range(cardinality)){
  		  Mtx[element][set] = 2.0 * ffd::random_distributions::Proba() -1.0;
      }
  	}

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

  	Rslt_minors = CalculateMinors(Mtx);
  	for (int minor : Range(1,Rslt.size())){
      for (int set : Range(cardinality)){
    		// std::cerr << minor << " " << set << " " << Rslt[minor][set] <<" " <<  Rslt_minors[minor][set] << std::endl;
  	 	  assert(std::abs(Rslt[minor][set] - Rslt_minors[minor][set]) <1e-10);
      }
    }
  }
}//namespace
