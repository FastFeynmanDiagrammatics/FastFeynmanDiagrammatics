namespace ffd::user_space::weighted_average::unit_test{

  void
  sum_two_means(){
    {
      ENumber<Real> x0{1., 1.};
      ENumber<Real> x1{2., 1.};
      WeightedMean<Real> w0(x0.eval());
      WeightedMean<Real> w1(x1.eval());
      w0 &= w1;
      auto const [val, err] = w0.eval();
      std::cerr << val << " " << err << "\n";
    }
    {
      ENumber<Real> x0{1., .1};
      ENumber<Real> x1{2., 1.};
      WeightedMean<Real> w0(x0.eval());
      WeightedMean<Real> w1(x1.eval());
      w0 &= w1;
      auto const [val, err] = w0.eval();
      std::cerr << val << " " << err << "\n";
    }
    
  }

}//namespace
