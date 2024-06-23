namespace ffd::nilpotent_polynomial::unit_test{

  void restrict_to_set(){
    
    using std::size;    using std::abs;
    
    
    const int order = 4;
    NilpotentPolynomial<long double> P(order);
    for(int j=0; j < 1<<order; ++j){
      P[j] = j+1;
    }

    
    auto Q0000 = P.restrict_to_set(0);

    assert(size(Q0000) == 1);
    static_assert(std::is_same<decltype(Q0000), decltype(P)>::value);
    assert(abs(Q0000[0]-P[0]) < abs(Q0000[0])*2*std::numeric_limits<long double>::epsilon());

    
    auto Q0001 = P.restrict_to_set(1);
    assert(size(Q0001) == 2);
    static_assert(std::is_same<decltype(Q0001), decltype(P)>::value);
    assert(abs(Q0001[0]-P[0]) < abs(Q0001[0])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q0001[1]-P[1]) < abs(Q0001[1])*2*std::numeric_limits<long double>::epsilon());

    
    auto Q0010 = P.restrict_to_set(0b0010);
    assert(size(Q0010) == 2);
    static_assert(std::is_same<decltype(Q0010), decltype(P)>::value);
    assert(abs(Q0010[0]-P[0]) < abs(Q0010[0])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q0010[1]-P[0b0010]) < abs(Q0010[1])*2*std::numeric_limits<long double>::epsilon());

    
    auto Q0100 = P.restrict_to_set(0b0100);
    assert(size(Q0100) == 2);
    static_assert(std::is_same<decltype(Q0100), decltype(P)>::value);
    assert(abs(Q0100[0]-P[0]) < abs(Q0100[0])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q0100[1]-P[0b0100]) < abs(Q0100[1])*2*std::numeric_limits<long double>::epsilon());


    auto Q1000 = P.restrict_to_set(0b1000);
    assert(size(Q1000) == 2);
    static_assert(std::is_same<decltype(Q1000), decltype(P)>::value);
    assert(abs(Q1000[0]-P[0]) < abs(Q1000[0])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1000[1]-P[0b1000]) < abs(Q1000[1])*2*std::numeric_limits<long double>::epsilon());


    auto Q0011 = P.restrict_to_set(0b0011);
    assert(size(Q0011) == 1<<2);
    static_assert(std::is_same<decltype(Q0011), decltype(P)>::value);
    assert(abs(Q0011[0]-P[0]) < abs(Q0011[0])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q0011[1]-P[0b0001]) < abs(Q0011[1])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q0011[0b10]-P[0b0010]) < abs(Q0011[0b10])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q0011[0b11]-P[0b0011]) < abs(Q0011[0b11])*2*std::numeric_limits<long double>::epsilon());

    
    auto Q0101 = P.restrict_to_set(0b0101);
    assert(size(Q0101) == 1<<2);
    static_assert(std::is_same<decltype(Q0101), decltype(P)>::value);
    assert(abs(Q0101[0]-P[0]) < abs(Q0101[0])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q0101[1]-P[0b0001]) < abs(Q0101[1])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q0101[0b10]-P[0b0100]) < abs(Q0101[0b10])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q0101[0b11]-P[0b0101]) < abs(Q0101[0b11])*2*std::numeric_limits<long double>::epsilon());


    auto Q1001 = P.restrict_to_set(0b1001);
    assert(size(Q1001) == 1<<2);
    static_assert(std::is_same<decltype(Q1001), decltype(P)>::value);
    assert(abs(Q1001[0]-P[0]) < abs(Q1001[0])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1001[1]-P[0b0001]) < abs(Q1001[1])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1001[0b10]-P[0b1000]) < abs(Q1001[0b10])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1001[0b11]-P[0b1001]) < abs(Q1001[0b11])*2*std::numeric_limits<long double>::epsilon());


    auto Q1010 = P.restrict_to_set(0b1010);
    assert(size(Q1010) == 1<<2);
    static_assert(std::is_same<decltype(Q1010), decltype(P)>::value);
    assert(abs(Q1010[0]-P[0]) < abs(Q1010[0])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1010[1]-P[0b0010]) < abs(Q1010[1])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1010[0b10]-P[0b1000]) < abs(Q1010[0b10])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1010[0b11]-P[0b1010]) < abs(Q1010[0b11])*2*std::numeric_limits<long double>::epsilon());


    auto Q1100 = P.restrict_to_set(0b1100);
    assert(size(Q1100) == 1<<2);
    static_assert(std::is_same<decltype(Q1100), decltype(P)>::value);
    assert(abs(Q1100[0]-P[0]) < abs(Q1100[0])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1100[1]-P[0b0100]) < abs(Q1100[1])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1100[0b10]-P[0b1000]) < abs(Q1100[0b10])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1100[0b11]-P[0b1100]) < abs(Q1100[0b11])*2*std::numeric_limits<long double>::epsilon());


    auto Q1110 = P.restrict_to_set(0b1110);
    assert(size(Q1110) == 1<<3);
    static_assert(std::is_same<decltype(Q1110), decltype(P)>::value);
    assert(abs(Q1110[0b000]-P[0b0000]) < abs(Q1110[0b000])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1110[0b001]-P[0b0010]) < abs(Q1110[0b001])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1110[0b010]-P[0b0100]) < abs(Q1110[0b010])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1110[0b011]-P[0b0110]) < abs(Q1110[0b011])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1110[0b100]-P[0b1000]) < abs(Q1110[0b100])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1110[0b101]-P[0b1010]) < abs(Q1110[0b101])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1110[0b110]-P[0b1100]) < abs(Q1110[0b110])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1110[0b111]-P[0b1110]) < abs(Q1110[0b111])*2*std::numeric_limits<long double>::epsilon());

    
    auto Q1101 = P.restrict_to_set(0b1101);
    assert(size(Q1101) == 1<<3);
    static_assert(std::is_same<decltype(Q1101), decltype(P)>::value);
    assert(abs(Q1101[0b000]-P[0b0000]) < abs(Q1101[0b000])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1101[0b001]-P[0b0001]) < abs(Q1101[0b001])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1101[0b010]-P[0b0100]) < abs(Q1101[0b010])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1101[0b011]-P[0b0101]) < abs(Q1101[0b011])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1101[0b100]-P[0b1000]) < abs(Q1101[0b100])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1101[0b101]-P[0b1001]) < abs(Q1101[0b101])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1101[0b110]-P[0b1100]) < abs(Q1101[0b110])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1101[0b111]-P[0b1101]) < abs(Q1101[0b111])*2*std::numeric_limits<long double>::epsilon());

    
    auto Q1011 = P.restrict_to_set(0b1011);
    assert(size(Q1011) == 1<<3);
    static_assert(std::is_same<decltype(Q1011), decltype(P)>::value);
    assert(abs(Q1011[0b000]-P[0b0000]) < abs(Q1011[0b000])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1011[0b001]-P[0b0001]) < abs(Q1011[0b001])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1011[0b010]-P[0b0010]) < abs(Q1011[0b010])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1011[0b011]-P[0b0011]) < abs(Q1011[0b011])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1011[0b100]-P[0b1000]) < abs(Q1011[0b100])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1011[0b101]-P[0b1001]) < abs(Q1011[0b101])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1011[0b110]-P[0b1010]) < abs(Q1011[0b110])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q1011[0b111]-P[0b1011]) < abs(Q1011[0b111])*2*std::numeric_limits<long double>::epsilon());


    auto Q0111 = P.restrict_to_set(0b0111);
    assert(size(Q0111) == 1<<3);
    static_assert(std::is_same<decltype(Q0111), decltype(P)>::value);
    assert(abs(Q0111[0b000]-P[0b0000]) < abs(Q0111[0b000])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q0111[0b001]-P[0b0001]) < abs(Q0111[0b001])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q0111[0b010]-P[0b0010]) < abs(Q0111[0b010])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q0111[0b011]-P[0b0011]) < abs(Q0111[0b011])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q0111[0b100]-P[0b0100]) < abs(Q0111[0b100])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q0111[0b101]-P[0b0101]) < abs(Q0111[0b101])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q0111[0b110]-P[0b0110]) < abs(Q0111[0b110])*2*std::numeric_limits<long double>::epsilon());
    assert(abs(Q0111[0b111]-P[0b0111]) < abs(Q0111[0b111])*2*std::numeric_limits<long double>::epsilon());


    auto Q1111 = P.restrict_to_set(0b1111);
    assert(size(Q1111) == 1<<4);
    static_assert(std::is_same<decltype(Q1111), decltype(P)>::value);
    for(int j=0; j < 1<<4; ++j){
      assert(abs(Q1111[j]-P[j]) < abs(Q1111[j])*2*std::numeric_limits<long double>::epsilon());
    }

    
  }

}//namespace
