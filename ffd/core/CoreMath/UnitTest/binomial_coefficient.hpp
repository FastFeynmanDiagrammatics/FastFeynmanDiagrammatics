namespace ffd::core_math::unit_test{

  void binomial_coefficient(){
    int constexpr n = 10;
    std::array<int, n+1> counter; counter.fill(0);
    for(BinaryInt S=0; S<(1<<n); ++S){
      counter[__builtin_popcount(S)]++;
    }

    // for(int k=0; k<=n; ++k){
    //   std::cerr<<counter[k]<<" "<<BinomialCoefficient(n, k)<<std::endl;
    // }

    {
      int constexpr k = 0;
        //     std::cerr<<counter[k]<<" "<<
	// BinomialCoefficient<n, k>()<<
	// std::endl;
      assert((counter[k] == BinomialCoefficient<n, k>() ));
      assert(( counter[k] == BinomialCoefficient(n, k) ));
    }

        {
      int constexpr k = 1;
      // std::cerr<<counter[k]<<" "<<
      // 	BinomialCoefficient<n, k>()<<
      // 	std::endl;
      assert((counter[k] == BinomialCoefficient<n, k>() ));
            assert(( counter[k] == BinomialCoefficient(n, k) ));
    }

	    {
      int constexpr k = 2;
      // std::cerr<<counter[k]<<" "<<
      // 	BinomialCoefficient<n, k>()<<
      // 	std::endl;
      assert((counter[k] == BinomialCoefficient<n, k>() ));
      assert(( counter[k] == BinomialCoefficient(n, k) ));
    }

	    {
      int constexpr k = 3;
// [[maybe_unused]]      std::array<int, BinomialCoefficient<n, k>()> X;
//       std::cerr<<counter[k]<<" "<<
// 	BinomialCoefficient<n, k>()<<
// 	std::endl;
      assert((counter[k] == BinomialCoefficient<n, k>() ));
            assert(( counter[k] == BinomialCoefficient(n, k) ));
    }

	    {
      int constexpr k = 4;
[[maybe_unused]]      std::array<int, BinomialCoefficient<n, k>()> X;
      // std::cerr<<counter[k]<<" "<<
      // 	BinomialCoefficient<n, k>()<<
      // 	std::endl;
      assert((counter[k] == BinomialCoefficient<n, k>() ));
            assert(( counter[k] == BinomialCoefficient(n, k) ));
    }

	    	    {
      int constexpr k = 5;
[[maybe_unused]]      std::array<int, BinomialCoefficient<n, k>()> X;
      // std::cerr<<counter[k]<<" "<<
      // 	BinomialCoefficient<n, k>()<<
      // 	std::endl;
      assert((counter[k] == BinomialCoefficient<n, k>() ));
            assert(( counter[k] == BinomialCoefficient(n, k) ));
    }

		    	    	    {
				      int constexpr k = 6;
[[maybe_unused]]      std::array<int, BinomialCoefficient<n, k>()> X;
      // std::cerr<<counter[k]<<" "<<
      // 	BinomialCoefficient<n, k>()<<
      // 	std::endl;
      assert((counter[k] == BinomialCoefficient<n, k>() ));
            assert(( counter[k] == BinomialCoefficient(n, k) ));
    }

				    		    	    	    {
				      int constexpr k = 7;
[[maybe_unused]]      std::array<int, BinomialCoefficient<n, k>()> X;
      // std::cerr<<counter[k]<<" "<<
      // 	BinomialCoefficient<n, k>()<<
      // 	std::endl;
      assert((counter[k] == BinomialCoefficient<n, k>() ));
                  assert(( counter[k] == BinomialCoefficient(n, k) ));
    }

		    	    	    {
				      int constexpr k = 8;
[[maybe_unused]]      std::array<int, BinomialCoefficient<n, k>()> X;
      // std::cerr<<counter[k]<<" "<<
      // 	BinomialCoefficient<n, k>()<<
      // 	std::endl;
      assert((counter[k] == BinomialCoefficient<n, k>() ));
                  assert(( counter[k] == BinomialCoefficient(n, k) ));
    }
		    	    	    {
				      int constexpr k = 9;
[[maybe_unused]]      std::array<int, BinomialCoefficient<n, k>()> X;
      // std::cerr<<counter[k]<<" "<<
      // 	BinomialCoefficient<n, k>()<<
      // 	std::endl;
      assert((counter[k] == BinomialCoefficient<n, k>() ));
                  assert(( counter[k] == BinomialCoefficient(n, k) ));
    }


				    		    	    	    {
				      int constexpr k = 10;
[[maybe_unused]]      std::array<int, BinomialCoefficient<n, k>()> X;
      // std::cerr<<counter[k]<<" "<<
      // 	BinomialCoefficient<n, k>()<<
      // 	std::endl;
      assert((counter[k] == BinomialCoefficient<n, k>() ));
                  assert(( counter[k] == BinomialCoefficient(n, k) ));
    }


  }

}//namespace
