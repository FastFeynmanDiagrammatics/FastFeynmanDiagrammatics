namespace ffd::bubbles::unit_test{

  template<bool IsLadder=true, bool EliminatingTadpoles=true>
  void gamma0_square_lattice(Real beta,
			     Real mu0,
			     Real U,
			     int L,
			     Real alpha=0.6,
			     Real precision=1e-4){



    using namespace std;

    using chebyshev_t = ffd::chebyshev_fft::ChebyshevFFT;

    using chebyshev_old_t = chebyshev_t;
      //ffd::chebyshev_polynomial::ChebyshevPolynomial<Real>;

    
    std::vector<std::array<chebyshev_t, 1>> G0(L*L);;
    std::vector<chebyshev_old_t> G0_OLD(L*L);
    for(int y=0; y<L; ++y){
      for(int x=0; x<L; ++x){
	auto G0_tau =
	  [mu0, L, beta, x, y](Real tau){
	    return G0_square_lattice(x,
				     y,
				     tau,
				     beta,
				     mu0,
				     L);
	  };

	
	G0[x+y*L][0] = chebyshev_t(G0_tau, {0, beta}, precision);

	
	G0_OLD[x+y*L] = chebyshev_old_t(G0_tau, {0, beta}, precision);
	
      }
    }


    int n_iter_GP = 30;
    ffd::user_space::Timer OLD_time;
    OLD_time.ini();
    ffd::rpa_ladder::RPA_Ladder<IsLadder, EliminatingTadpoles,
				chebyshev_old_t> bubbles_OLD(G0_OLD, beta, U, alpha, precision);
    // bubbles_OLD.IterateUntilRequestedPrecision();
    using ffd::vector_range::Range;
    for( int j: Range(n_iter_GP) ){
      bubbles_OLD.IterateG();
      bubbles_OLD.IterateP();
    }
    std::cerr<<"OLD time = "<<OLD_time.elapsed()<<"s\n";
    

    
    ffd::imaginary_time::PeriodicImaginaryTime BetaStrip(beta);
    ffd::lattice::HypercubicLattice<2> Square({L, L});
    // ffd::lattice::HypercubicLattice<1> Square({L});
    

    RPA_scheme scheme;
    if(IsLadder){
      scheme = RPA_scheme::pp;
    }else{
      scheme = RPA_scheme::ph;
    }
    tadpoles tad;
    if(EliminatingTadpoles){
      tad = tadpoles::no;
    }else{
      tad = tadpoles::yes;
    }

    

    ffd::user_space::Timer NEW_time;
    NEW_time.ini();
    Bubbles bubbles(BetaStrip, Square, {G0, G0}, U, scheme, tad);
    bubbles.AbsolutePrecision = precision;
    bubbles.damping_alpha = alpha;
    // bubbles.IterateAll();
    for( int j: Range(n_iter_GP) ){
      bubbles.IterateG();
      bubbles.IterateP();
    }
    std::cerr<<"NEW VERSION time = "<<NEW_time.elapsed()<<"s\n";


    int const N_samples_tau = 10;
    for(int x=0; x<L; ++x){
      for(int y=0; y<L; ++y){
	for(int k=0; k<N_samples_tau; ++k){
	  Real tau = (k+.3)*beta/N_samples_tau;
	  //ffd::user_space::Proba()*beta;
	  Real G_v = bubbles.G[0][x+L*y][0](tau);
	  Real P_v = bubbles.P[0][x+L*y](tau);
	  std::cerr<<x<<" "<<y<<": "<<G_v<<" "<<G_v-bubbles_OLD.G[x+L*y](tau)<<std::endl;
	  std::cerr<<x<<" "<<y<<": "<<P_v<<" "<<P_v-bubbles_OLD.P[x+L*y](tau)<<std::endl;
	}
      }
    }
    
  }

}//namespace
