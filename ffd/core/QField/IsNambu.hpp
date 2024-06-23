namespace ffd::user_space{

  bool
  IsNambu(QField x){
    using namespace q_field;
    return value_bit(x, nambu_bit);
  }

}//namespace
