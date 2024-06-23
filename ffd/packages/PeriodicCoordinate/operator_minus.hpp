namespace ffd::user_space{

	template<
   	   int d,
   	  int rec,
   	   typename... args>
    void
     operator_minus_rec(ffd::periodic_coordinate::PeriodicCoordinates<d, args...>& C1, std::tuple<args...> shift){
       if constexpr(rec < sizeof...(args)){
   		component<rec>(C1) -= std::get<rec>(shift);
   		operator_minus_rec<d,rec+1,args...>(C1,shift);
   	}
     }

     template<
    	   int d,
    	   typename... args>
     auto
      operator-(ffd::periodic_coordinate::
		PeriodicCoordinates<d, args...> const& C1,
		std::tuple<args...> shift){
   	   auto C2 = C1;
   	   operator_minus_rec<d,0,args...>(C2,shift);
   	   return C2;
      }

}//namespace
