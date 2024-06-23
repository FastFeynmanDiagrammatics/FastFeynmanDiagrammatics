namespace ffd::ranked_zeta::unit_test{

  void test_extend(){
    using NilpotentPolynomial = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Full>;
    using EvenNilPoly = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Even>;
    using namespace user_space;
    using std::size;

    const int order = 10;
    int set_1=562, set_2=687, set_3=686;

    NilpotentPolynomial P1(4), P2(7);
    for(int j=0; j < 1<<4; ++j){
      P1[j] = std::sin(0.7+4.4*j);
    }
    for(int j=0; j < 1<<7; ++j){
      P2[j] = std::sin(1.7+2.1*j);
    }

    RankedZeta Z1 = P1;
    RankedZeta Z2 = P2;

    P1=P1.extend_to_set(set_1, order);
    P2=P2.extend_to_set(set_2, order);
    P1=P1.shift((1<<order)-1-set_1, order);
    P2=P2.shift((1<<order)-1-set_2, order);

    auto P3 = P1+P2;

    Z1=Z1.extend_to_set(set_1, order);
    Z2=Z2.extend_to_set(set_2, order);

    Z1=Z1.shift((1<<order)-1-set_1, order);
    Z2=Z2.shift((1<<order)-1-set_2, order);

    auto Z3 = Z1+Z2;

    NilpotentPolynomial P3_new = Z3;

    for(BinaryInt j=0; j < P3.size(); ++j){
      // std::cerr << j << " " <<  P3[j] << " " << P3_new[j] << " " << std::abs(P3[j] - P3_new[j]) << std::endl;
      if (P3[j]!=0. || P3_new[j]!=0.) assert(std::abs(P3[j] - P3_new[j]) <= std::max(P3[j] * 1.e-10, 1.e-10) );
    }

    EvenNilPoly P1_e(4), P2_e(6);
    for(int j=0; j < 1<<4; ++j){
      if (__builtin_parity(j)==0) P1_e[j] = std::sin(0.7+4.4*j);
    }
    for(int j=0; j < 1<<6; ++j){
      if (__builtin_parity(j)==0) P2_e[j] = std::sin(0.7+4.4*j);
    }

    // std::cerr << "here" << std::endl;

    EvenZeta
      Z1_e = P1_e;
    EvenZeta Z2_e = P2_e;

    // std::cerr << "here2" << std::endl;

    auto P1_ext=P1_e.extend_to_set(set_1, order);
    auto P2_ext =P2_e.extend_to_set(set_3, order);
    P1_ext=P1_ext.shift((1<<order)-1-set_1, order);
    P2_ext=P2_ext.shift((1<<order)-1-set_3, order);

    auto P3_ext = P1_ext+P2_ext;

    // std::cerr << "created" << std::endl;
    // 		std::cerr << Z1_e.rank_size << " " << Z2_e.rank_size<< " " << Z1_e.cardinality<< " " << Z2_e.cardinality << std::endl;
    //      std::cerr << Z1_e.size() << " " << Z2_e.size()<< std::endl;

    Z1_e=Z1_e.extend_to_set(set_1, order);
    Z2_e=Z2_e.extend_to_set(set_3, order);

    // std::cerr << "extended" << std::endl;
    // 		std::cerr << Z1_e.rank_size << " " << Z2_e.rank_size<< " " << Z1_e.cardinality<< " " << Z2_e.cardinality << std::endl;
    //          std::cerr << Z1_e.size() << " " << Z2_e.size()<< std::endl;

    Z1_e=Z1_e.shift((1<<order)-1-set_1, order);
    Z2_e=Z2_e.shift((1<<order)-1-set_3, order);
    // std::cerr << "shifted " << std::endl;
    // std::cerr << Z1_e.rank_size << " " << Z2_e.rank_size<< " " << Z1_e.cardinality<< " " << Z2_e.cardinality << std::endl;
    //  std::cerr << Z1_e.size() << " " << Z2_e.size()<< std::endl;
    auto Z3_e = Z1_e+Z2_e;

    // std::cerr << "combined" << std::endl;

    EvenNilPoly P3_e_new = Z3_e;

    // std::cerr << "here7" << std::endl;

    for(BinaryInt j=0; j < P3_ext.size(); ++j){
      // std::cerr << j << " " <<  P3[j] << " " << P3_new[j] << " " << std::abs(P3[j] - P3_new[j]) << std::endl;
      if (P3_ext[j]!=0. || P3_e_new[j]!=0.) assert(std::abs(P3_ext[j] - P3_e_new[j]) <=
						   std::max(P3_ext[j] * 1.e-10, 1.e-10) );
    }
  }

}//namespace
