namespace ffd::user_space::imaginary_time_lattice_propagator::unit_test{

  template<int d>
  
  auto
  assertable_Nambu_propagator(std::array<int, d> L,
			      std::array<Real, 2> Mu,
			      std::array<std::array<Real, d>, 2> t_spin_j,
			      Real Delta,
			      Real beta){
    bool IsOk = true;


    using ffd::vector_range::Range;


    
    ffd::lattice::HypercubicLattice<d> Hypercubic(L);
    ffd::imaginary_time::PeriodicImaginaryTime BetaStrip(beta);

    

    auto O = CreateCoordinates(BetaStrip,
			       Hypercubic);
    component<0>(O) = 1e-15;

    
    QuadraticAction<Real> S0;
    for( int spin: {-1, 1} ){
      S0 += Bar(Psi_(spin))(O)*Psi_(spin)(O) * ( - Mu[(1+spin)/2] );
    }


    auto hopping = hypercubic_hopping<d, 0>(O, t_spin_j);
    S0 += hopping;



    auto anomalous = Psi_(-1)(O)*Psi_(1)(O)*(-Delta);
    S0 += anomalous + HermitianConjugate(anomalous);

    auto S0_Nambu = Nambu(S0);

    
    auto G0 = ImaginaryTimeLatticePropagator(BetaStrip,
					     Hypercubic,
					     S0);

    
    // G0.BuildChebyshev();
    
    
    auto G0_Nambu = ImaginaryTimeLatticePropagator(BetaStrip,
						   Hypercubic,
						   S0_Nambu);
    
    
    NambuPropagatorExact<d> G0_exact;

    
    G0_exact.L = L;
    G0_exact.Mu = Mu;
    G0_exact.t_spin_j = t_spin_j;
    G0_exact.beta = beta;
    G0_exact.Delta = Delta;


    
    std::array<ffd::quantum_field::QuantumField, 2>
      psi_array;
    psi_array[0] = Psi_(1)[0];
    psi_array[1] = Bar(Psi_(-1)[0]);

    
    int N_samples_tau = 4;
    for( int x: Range(3) ){
      for( int t: Range(N_samples_tau) ){
	Real tau = (t+.2)*beta/N_samples_tau;

	
	for( int index0: {0, 1} ){
	  for( int index1: {0, 1} ){
	    auto psi = psi_array[index0];
	    auto X = O;
	    component<0>(X) += tau;
	    component<1>(X) += x;
	    auto pair_psi  = std::pair{psi, std::any(X)};
	    auto pair_psi_n = std::pair{Nambu(psi), std::any(X)};
	    
	    auto bpsi = Bar(psi_array[index1]);
	    auto pair_bpsi = std::pair{bpsi, std::any(O)};
	    auto pair_bpsi_n = std::pair{Nambu(bpsi), std::any(O)};

	    std::array<int, d> R;
	    R.fill(0);
	    R[0] = x;

	    Real G0_value = G0({{pair_psi, pair_bpsi}});
	    // std::cerr<<"("<<index0<<", "<<index1<<") "<<G0_value<<" "<<
	    //   G0_Nambu({{pair_psi_n, pair_bpsi_n}})-G0_value<<" "<<
	    //   -(2*index0-1)*(2*index1-1)*G0_exact( tau, R, {{index0, index1}} )-G0_value<<std::endl;
	    IsOk = IsOk &&
	      std::abs(G0_Nambu({{pair_psi_n, pair_bpsi_n}}) -
		       G0_value) <
	      20*std::numeric_limits<Real>::epsilon();

	    IsOk = IsOk &&
	      std::abs(-(2*index0-1)*(2*index1-1)*
		       G0_exact( tau, R, {{index0, index1}} ) -
		       G0_value) <
	      1e-4;
	  }
	}
      }
    }    
    // std::cerr<<"\n\n";
    
    return IsOk;
  }


}//namespace
