namespace ffd::determinant{

  template<typename matrix_t>

  typename std::decay<decltype(std::declval<matrix_t>()[0])>::type
  Determinant(matrix_t&& A);
	
	
}//namespace
