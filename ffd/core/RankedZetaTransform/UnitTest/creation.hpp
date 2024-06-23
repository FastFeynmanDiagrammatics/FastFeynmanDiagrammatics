namespace ffd::ranked_zeta::unit_test{

	void creation(){
      RankedZeta<Real> Z;

      assert(!Z.NotZero);
      assert(Z.size() == 1);
      assert(std::abs(Z[0]) < 10*std::numeric_limits<Real>::min() );

      const int order = 5;
      RankedZeta<float> Z2(order);

      assert(!Z2.NotZero);
      Z2[0] = 0.;
      assert(Z2.NotZero);

      assert(Z2.size() == (1<<order)*(order+1));
      assert(Z2.get_cardinality() == (1<<order));
      assert(Z2.get_rank_size() == order);

      for(int j=0; j < 1<<order; ++j){
        assert(std::abs(Z2[j]) < 10*std::numeric_limits<Real>::min() );
      }

      const long double z = 3.12312312312392438232323423748237l;
      RankedZeta<Real> Z3 = z;
      assert(Z3.NotZero);
      assert(Z3.size() == 1);
      assert(std::abs(Z3[0]-z) < 10*std::numeric_limits<Real>::epsilon() );

      std::vector<RankedZeta<Real>> Z_vec;
    }


}//namespace
