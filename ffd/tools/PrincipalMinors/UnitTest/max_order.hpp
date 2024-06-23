namespace ffd::principal_minors::unit_test{

  void test_max_order(){
    using ffd::vector_range::Range;

    std::vector<Real> Mtx(6*6,0.0);
    std::vector<Real> Rslt_minors((1<<6),0.0);
    std::vector<Real> Rslt_minors_step((1<<3),0.0);

    for (int element : Range(Mtx.size())){
      Mtx[element] = cos(element);
    }

    const uint max_ord = 2;

    Rslt_minors = DeterminantPrincipalMinors<max_ord>(Mtx);
    for (int minor : Range(Rslt_minors.size())){
      // std::cerr << minor << " " << Rslt_minors[minor] << std::endl;
      if (__builtin_popcount(minor)>max_ord){
        assert(Rslt_minors[minor]==0);
      }
    }

    Rslt_minors_step = DeterminantPrincipalMinors<max_ord>(Mtx,2);
    for (int minor : Range(Rslt_minors_step.size())){
      if (__builtin_popcount(minor)>max_ord){
        assert(Rslt_minors_step[minor]==0);
      }
    }
  }

}//namespace
