namespace ffd::lattice{

  template<typename LatticeType>
  
  std::vector<int>
  
  GetLinearSizes(LatticeType const& Lattice1){
    std::vector<int> L(LatticeType::dimension, 0);

    
    auto const x_linear_sizes = CreateCoordinates(Lattice1);
    ffd::lattice::
      get_linear_sizes<LatticeType::dimension,
		       0>(x_linear_sizes, L);

    
    return L;
  }


 
}//namespace
