namespace ffd::chebyshev_fft{

  ChebyshevFFT ChebyshevFFT::Derivative(int order_derivative) const{
    assert(( order_derivative >= 0 ));
    if(order_derivative == 0){
      return *this;
    }else{
      return (this->Derivative()).Derivative(order_derivative-1);
    }
  }


}//namespace
