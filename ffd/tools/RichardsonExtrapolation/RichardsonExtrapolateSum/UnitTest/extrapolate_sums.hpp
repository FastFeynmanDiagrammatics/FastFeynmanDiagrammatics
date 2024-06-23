

namespace ffd::richardson_extrapolation::sums::unit_test{

  void extrapolate_sums(){
    using namespace std;
    Real epsilon = 1e-10;
    auto value = RichardsonExtrapolateSum([](int n)->Real{return 1./(1+n*n);},
					  {"-Infinity", "Infinity"}, epsilon);
  }


}//namespace 
