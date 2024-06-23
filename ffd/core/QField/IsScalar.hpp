namespace ffd::user_space{

  bool
  IsScalar(QField x){
    using namespace q_field;
    return !value_bit(x, scalar_bit);
  }

}//namespace
