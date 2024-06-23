

namespace ffd::rpa_ladder::unit_test{

  void UnitTest(){
    Real beta = 1.1413241, mu0 = -.3, U = 8, alpha = .6, precision = 1e-8;

    
    bool const Ladder = true;
    bool const EliminatingTadpoles = true;


    
    assert((
	    gamma0_atom<Ladder,
	    !EliminatingTadpoles>(
				  beta,
				  mu0,
				  U,
				  alpha,
				  precision
				  )
	    ));


    

    assert((
	    sunset_atom<Ladder,
	    !EliminatingTadpoles>(
				  beta,
				  mu0,
				  U,
				  alpha,
				  precision)
	    ));

    
    

    int L = 4;
    gamma0_square_lattice<Ladder,
			  !EliminatingTadpoles>(
						beta,
						mu0,
						U=4,
						L,
						alpha=.3,
						precision=1e-5);


    
    gamma0_square_lattice<Ladder,
			  EliminatingTadpoles>(
						beta,
						mu0,
						U=4,
						L,
						alpha=.3,
						precision=1e-5);


    

    
    assert((
	    check_if_atomic_lobster_is_cooked<Ladder,
	    !EliminatingTadpoles>(beta=1.12,
				  mu0=-.3,
				  U=8,
				  alpha=.6,
				  precision=1e-4)
	    ));


    
    // assert((
    // 	    check_if_atomic_lobster_is_cooked<Ladder,
    // 	    !EliminatingTadpoles>(beta=2.12,
    // 				  mu0=-.6,
    // 				  U=4,
    // 				  alpha=.6,
    // 				  precision=1e-12)
    // 	    ));


    
    // assert((
    // 	    check_if_atomic_lobster_is_cooked<Ladder,
    // 	    !EliminatingTadpoles>(beta=3.12,
    // 				  mu0=1.2,
    // 				  U=10,
    // 				  alpha=.6,
    // 				  precision=1e-12)
    // 	    ));

	


    check_if_lobster_is_cooked<Ladder,
			       EliminatingTadpoles>(beta=1,
						    mu0=-.3,
						    U=4,
						    L=2,
						    alpha=.6,
						    precision=1e-4);

    
    
  }
  
  
}//namespace

