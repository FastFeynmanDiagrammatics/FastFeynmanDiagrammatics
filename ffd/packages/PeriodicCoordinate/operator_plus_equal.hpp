
namespace ffd::user_space{

     template<
    	   int d,
    	   typename... args>
     auto &
      operator+=(ffd::periodic_coordinate::PeriodicCoordinates<d, args...> & C1, std::tuple<args...> shift){
   	   operator_plus_rec<d,0,args...>(C1,shift);
   	   return C1;
      }

}//namespace
