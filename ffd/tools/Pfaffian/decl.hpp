namespace ffd::gauss_pfaffian{

  template<typename antisymmetric_matrix_t>

  typename std::decay<decltype(std::declval<antisymmetric_matrix_t>()[0])>::type
  Pfaffian(antisymmetric_matrix_t&& A);

}//namespace
