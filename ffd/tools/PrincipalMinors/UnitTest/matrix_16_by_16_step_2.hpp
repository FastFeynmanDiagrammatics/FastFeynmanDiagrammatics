namespace ffd::principal_minors::unit_test{

  void matrix_16_by_16_step_2(){
    using ffd::vector_range::Range;
    using ffd::set_theory::BinarySpread;

    std::vector<Real> Mtx(16*16,0.0);
    // std::vector<Real> Rslt_minors((1<<16),0.0);
    std::vector<Real> Rslt_dets((1<<8),0.0);
    std::vector<Real> Rslt_minors_step((1<<8),0.0);
    std::vector<Real> Rslt_minors_2n((1<<8),0.0);
    // for (int cnt=0; cnt<(1<<1); cnt++)
    {
      for (int element : Range(Mtx.size())){
        // Mtx[element] = 2.0 * ffd::random_distributions::Proba() -1.0;
    		 Mtx[element] = 2.0 * cos(element+2.31542) + abs(sin(element*2.3193)) -0.152/(double)(1+element);
      }

//      ffd::math_tools::print_matrix(Mtx);

      // Rslt_minors = PrincipalMinors(Mtx);
      Rslt_dets = DeterminantPrincipalMinors(Mtx,2);
      Rslt_minors_step = PrincipalMinors(Mtx,2);
      // Rslt_minors_old = PMOLD(Mtx,2);
      Rslt_minors_2n = PrincipalMinors2n(Mtx,2);
      for (int minor : Range(Rslt_minors_step.size())){
        // if (std::abs(Rslt_minors_step[minor] - Rslt_dets[minor]) > abs(Rslt_dets[minor])/1.e-10 ||
        //     std::abs(Rslt_minors_2n[minor]   - Rslt_dets[minor]) > abs(Rslt_dets[minor])/1.e-10)
        // std::cerr << minor << " " << Rslt_minors_step[minor]
        //                    << " " << Rslt_minors_2n[minor]
        //             			 << " " << Rslt_dets[minor]
        //             			 << " " << std::abs(Rslt_minors_step[minor] - Rslt_dets[minor])
        //             			 << " " << std::abs(Rslt_minors_2n[minor] - Rslt_dets[minor])
        //             			 << " " << std::endl;
        assert(std::abs(Rslt_minors_step[minor] - Rslt_dets[minor]) < abs(Rslt_dets[minor])/1.e-10);
        assert(std::abs(Rslt_minors_2n[minor] - Rslt_dets[minor]) < abs(Rslt_dets[minor])/1.e-10);
      }
    }
  }

}//namespace
