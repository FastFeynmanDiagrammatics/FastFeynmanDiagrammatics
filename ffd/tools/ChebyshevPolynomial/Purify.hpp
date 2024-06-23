namespace ffd::chebyshev_polynomial{

  template<typename Field>

  void
  ChebyshevPolynomial<Field>::Purify(Real precision){

    
    int new_size = size(*this);
    for( ;
	 new_size > 0 && std::abs(this->Coef[new_size-1]) < precision ;
	 --new_size ){}

    
    this->Coef.resize(new_size);
    this->Coef.shrink_to_fit();
  }

}//namespace
