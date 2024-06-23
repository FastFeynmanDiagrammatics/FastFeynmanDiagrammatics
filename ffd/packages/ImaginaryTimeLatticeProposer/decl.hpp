namespace ffd::imaginary_time_lattice_proposer{

  template<typename imaginary_time_t,
	   typename lattice_t>
  struct Proposer{
    imaginary_time_t imaginary_time_v;
    lattice_t lattice_v;


    ffd::chebyshev_fft::ChebyshevFFT time_function;
    ffd::chebyshev_fft::ChebyshevFFT time_cumulative;


    std::array<std::vector<Real>, lattice_t::number_atoms_unit_cell> space_function;
    std::array<std::vector<Real>, lattice_t::number_atoms_unit_cell> space_cumulative;
    using array_space_components_t =
      std::array<int, lattice_t::dimension+(lattice_t::number_atoms_unit_cell>1)>;
    std::vector<array_space_components_t>
    vector_lattice_sites_components;
    

    
    using coordinates_t =  decltype(ffd::user_space::CreateCoordinates(std::declval<imaginary_time_t>(),
								       std::declval<lattice_t>()));




    template<typename time_seed_function_t,
	     typename space_seed_function_t>
    Proposer(imaginary_time_t imaginary_time_v_,
	     lattice_t lattice_v_,
	     time_seed_function_t time_seed,
	     space_seed_function_t space_seed):
      imaginary_time_v(imaginary_time_v_),
      lattice_v(lattice_v_){
      create_time_cumulative(time_seed);
      create_space_cumulative(space_seed);
    }



    template<typename time_seed_function_t>
    void create_time_cumulative(time_seed_function_t);
    template<typename space_seed_function_t>
    void create_space_cumulative(space_seed_function_t);

    
    Real                     propose_time() const;
    array_space_components_t propose_space(int atom = 0) const;




    coordinates_t  propose(coordinates_t) const;
    coordinates_t  operator()(coordinates_t) const;//alias for propose( )
    

    
    Real operator()(coordinates_t, coordinates_t) const;


  };
	
	
}//namespace
