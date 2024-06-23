namespace ffd::bubbles{

  enum class RPA_scheme{ph, pp};

  enum class tadpoles{yes, no};

  
  template<typename ImaginaryTime_t,
	   typename Lattice_t>
  class Bubbles{
  public:
    using chebyshev_t = ffd::chebyshev_fft::ChebyshevFFT;
    using chebyshev_2d_t = ffd::chebyshev_2d_fft::Chebyshev2DFFT;
    
    static constexpr int d = Lattice_t::dimension;
    static constexpr int cell_size = Lattice_t::number_atoms_unit_cell;

    using Gamma_t = std::vector<std::array<chebyshev_t,
					   Lattice_t::number_atoms_unit_cell*
					   Lattice_t::number_atoms_unit_cell>>;
    using GF_t = std::array<Gamma_t, 2>;

    using Lobster_bravais_t = std::array<chebyshev_2d_t,
					 Lattice_t::number_atoms_unit_cell*
					 Lattice_t::number_atoms_unit_cell*
					 Lattice_t::number_atoms_unit_cell*
					 Lattice_t::number_atoms_unit_cell>;
    
    using Lobster_t = std::vector<Lobster_bravais_t>;
    
    ImaginaryTime_t ImagTime;
    Lattice_t Lattice;

    GF_t G0;
    GF_t G;

    Real beta;
    Real U;
    RPA_scheme const RPA;
    tadpoles tadpoles_choice;
    Real AbsolutePrecision = 1e-10;
    Real damping_alpha = 0.6;

    Gamma_t P;

    std::array<int, d> L_array;
    int L_product;
    
    Lobster_t Lobster;
    
    // (function/class/member) Lobster

    //   (function/class) Shrimp




    Bubbles(ImaginaryTime_t ImagTime_,
	    Lattice_t Lattice_,
	    GF_t G0_,
	    Real U_,
	    RPA_scheme RPA_ = RPA_scheme::pp,
	    tadpoles tadpoles_choice_ = tadpoles::no):
      ImagTime(ImagTime_),
      Lattice(Lattice_),
      G0(G0_),
      G(G0_),
      beta( Beta(ImagTime_) ),
      U(U_),
      RPA(RPA_),
      tadpoles_choice(tadpoles_choice_)
    {
      auto L = ffd::lattice::
	GetLinearSizes(Lattice);
      L_product = 1;
      for( int j: ffd::vector_range::Range(d) ){
	L_product *= L[j];
	L_array[j] = L[j];
      }


      P.resize(L_product);
      for( auto& p_el: P ){
	for( auto& p_el_el: p_el ){
	  p_el_el = chebyshev_t();
	}
      }
      
    }      


    using error_t = Real;
    

    error_t
    IterateP();
    error_t
    IterateG();
    void
    IterateAll();

    
    void
    CookLobster(Real precision = -1.);
    

  };
	
}//namespace


