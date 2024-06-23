namespace ffd::user_space{
  
  template<typename imaginarytime_t,
	   typename lattice_t,
	   typename actionfield_t,
	   typename propagatorfield_t = Real>
  
  struct ImaginaryTimeLatticePropagator{
    
    using CoordinatesType = decltype( CreateCoordinates(std::declval<imaginarytime_t>(),
							std::declval<lattice_t>() ) );
    
    static constexpr int d                      = lattice_t::dimension;
    static constexpr int number_atoms_unit_cell = lattice_t::number_atoms_unit_cell;
    
    imaginarytime_t ImagTimeOfModel;
    
    lattice_t LatticeOfModel;
    
    QuadraticAction<actionfield_t> Action;
    
    Real AbsolutePrecision;


    ffd::flat_map::FlatMap<
      std::tuple<std::array<ffd::quantum_field
			    ::QuantumField, 2>,
		 std::array<int, 2>,
		 std::array<int, d>>,
		 ffd::chebyshev_polynomial::
		 ChebyshevPolynomial<propagatorfield_t>>
    ChebyshevPropagator;
    

    bool UseChebyshevPolynomials = false;
    
    
    ffd::wick_block_oracle::WickBlockOracle BlockOracle;

    std::vector<std::optional<Real>>
    TwistedBoundaryConditionsVector;
    
    std::vector<Real> InverseBravaisMatrix;

    
    using BlockTerm =
      std::tuple<		      
      actionfield_t,
      std::array<std::array<Real, d>, 2>,
      std::array<int, 2>,
      std::array<bool, 2>>;
    
    
    std::vector<std::vector<BlockTerm>> BlocksTerms;
    
    
    std::vector<bool> IsAFermionicBlock;
    std::vector<int> BlockSize;
    std::vector<std::optional<bool>> NotAnomalousBlock;
    
    
    std::vector<
      std::vector<
	std::pair<std::vector<Real>,
		  std::vector<std::vector<Complex>>>>>
    EigenvaluesVectorsOfBlock1Momentum2;
    
    
    
    ImaginaryTimeLatticePropagator(imaginarytime_t imag_time_,
				   lattice_t lattice_,
				   QuadraticAction<actionfield_t> action_,
				   propagatorfield_t precision_and_propagator_type_ = 1e-10):
      ImagTimeOfModel(imag_time_),
      LatticeOfModel(lattice_),
      Action(action_),
      AbsolutePrecision(std::abs(precision_and_propagator_type_)),
      BlockOracle(action_),
      TwistedBoundaryConditionsVector(lattice_t::dimension)
    {
      ComputeInverseBravais();
      ComputeBlockTerms();
      ComputeEigenvaluesEigenvectors();
    }
    
    

    void ComputeInverseBravais();

  private:
    auto GetRealSpacePositions(std::any SpacetimePosition) const;
    auto NumberParticlesInBlock(int block_number) const;
    int GetParticleIndex(ffd::quantum_field::QuantumField) const;
    int GetAtomIndex(std::any SpacetimePosition) const;
    int CreateAtomParticleIndex(int atom_index, int particle_index) const;
    
  public:
    auto
    ReturnBlockTerm(int block_number) const;

  public:
    void ComputeBlockTerms();

  private:
    std::vector<int> GetLinearSizes() const;
    Real GetInverseTemperatureOfModel() const;

  public:    
    std::vector<std::vector<Complex>>
    ReturnBlockHamiltonian_k(std::vector<BlockTerm> const&
			     terms_quadratic_action_block,
			     bool is_fermionic,
			     int block_size,
			     std::array<Real, d>
			     wave_vector) const;

  public:
    

    void ComputeEigenvaluesEigenvectors();
    

    propagatorfield_t
    operator()(ffd::feynman_edge::
	       QuantumFieldPositions) const;

    
    void BuildChebyshev();

    
  };
  


}//namespace
