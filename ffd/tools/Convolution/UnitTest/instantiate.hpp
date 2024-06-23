namespace ffd::convolution::unit_test{

  void instantiate(){

    auto f= [](Real x){return x;};
    auto g= [](Real x){return x*x;};


    Convolution(f, g, {-1., 2.});

  }
  


}//namespace
