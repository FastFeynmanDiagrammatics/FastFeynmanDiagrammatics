namespace ffd::eigenvalues_vectors{

  template<typename Field>
  std::pair<std::vector<Real>, std::vector<std::vector<Field>>>
  EigenvaluesVectors(std::vector<std::vector<Field>> const& MatrixToBeDiagonalized_,
		     Real precision = -1.);
  	
	
}//namespace
