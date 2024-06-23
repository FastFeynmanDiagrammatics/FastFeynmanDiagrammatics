

namespace ffd::eigenvalues_vectors{

  auto PhaseMatrix(Real phase,
		   std::array<int, 2> components_to_be_dephased,
		   int size_matrix){
    auto phase_matrix =
      UnitMatrix<Complex>(size_matrix);

    
    auto const& k = components_to_be_dephased;
    assert( k[0] != k[1] );

    

    using ffd::user_space::expI;



    for(int j: {0, 1}){
      phase_matrix[ k[j] ][ k[j] ] = expI( (.5-j)*phase );
    }
    return phase_matrix;
  }
  

}//namespace
