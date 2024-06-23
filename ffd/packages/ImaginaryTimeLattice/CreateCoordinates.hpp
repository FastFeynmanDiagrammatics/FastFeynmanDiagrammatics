namespace ffd::user_space{


  template<int d, auto... args,
	   template<int, auto...> class lattice,
	   bool NotAntiperiodic = false>
  auto CreateCoordinates(ffd::imaginary_time::ImaginaryTime<NotAntiperiodic> const& Imag_time_,
			 lattice<d, args...> const& Lattice_){
    return ffd::merge(
		      CreateCoordinates<d>(Imag_time_),
		      CreateCoordinates(Lattice_)
		      );
  }

  
  template<int d, auto... args,
	   template<int, auto...> class lattice,
	   typename... tuple_t,
	   bool NotAntiperiodic = false>
  auto CreateCoordinates(ffd::imaginary_time::ImaginaryTime<NotAntiperiodic> const& Imag_Time_,
			 lattice<d, args...> const& Lattice_, std::tuple<tuple_t...> X){
    auto Y = CreateCoordinates(Imag_Time_,Lattice_);
    return Y+X;
  }


  
  template<int d, auto... args,
	 template<int, auto...> class lattice,
	 bool NotAntiperiodic = false>
  auto CreateCoordinates(ffd::imaginary_time::ImaginaryTime<NotAntiperiodic> const& Imag_Time_,
			 lattice<d, args...> const& Lattice_, std::size_t l){

    std::vector<decltype(CreateCoordinates(Imag_Time_, Lattice_))> Y(l);
    for (size_t i=0; i<l; ++l){
      Y[i] = CreateCoordinates(Imag_Time_,Lattice_);
    }
    return Y;
  }

}//namespace
