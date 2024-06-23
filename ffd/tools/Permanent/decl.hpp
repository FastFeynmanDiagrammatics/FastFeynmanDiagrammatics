namespace ffd::permanent{

  template<typename matrix_t>

  typename std::decay<decltype(std::declval<matrix_t>()[0])>::type
  Permanent(matrix_t&& A);
	
	
}//namespace
