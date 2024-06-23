namespace ffd::distributions{

 auto pade_gauss(Real dev = 1.){
     using std::abs;
     using std::exp;
     using namespace ffd::type_traits;
     
     auto f = [dev] (auto x){
     if constexpr (!is_container<decltype(x)>::value){
     auto y = abs(x)/dev;
     auto y2 = y*y;
     return 1./(1.+.5*y2*(1.+0.25*y2));
     }
    else{
     auto y = abs(x[0]);
     auto x2 = y*y;
     for (std::size_t j=1; j<size(x); ++j){
     y = abs(x[j]);
     x2 += y*y;
     }
     x2/=(dev*dev);
     return 1./(1.+.5*x2*(1.+0.25*x2));
    }
   };
   return f; 
 } 

}//namespace
