namespace ffd::distributions{

  auto gaussian(Real dev = 1.){
    using std::abs;
    using std::exp;
    using namespace ffd::type_traits;
     

    auto f =
      [dev] (auto x){
	if constexpr (!is_container<decltype(x)>::value){
	    auto y = abs(x);
	    return exp(-0.5*y*y/(dev*dev));
	  }
	else{
	  auto y = abs(x[0]);
	  auto x2 = y*y;
	  for (std::size_t j=1; j<size(x); ++j){
	    y = abs(x[j]);
	    x2 += y*y;
	  }
	  return exp(-0.5*x2/(dev*dev));       
	}
      };
    return f;
  } 

}//namespace
